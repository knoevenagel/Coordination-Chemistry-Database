# app.py
import csv
import os
import time
import sqlite3
from typing import Dict, List, Any, Optional

import ujson
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

# ======================
# 配置：5 个 CSV 与索引列（等值查询）
# ======================
DATASETS_CONFIG = [
    {
        "name": "l1",
        "csv_path": "./tmp/complex_data.csv",
        "index_cols": ["did"],
    },
    {
        "name": "l2",
        "csv_path": "./tmp/ligand_data.csv",
        "index_cols": ["ligand_did"],
    },
    {
        "name": "l3",
        "csv_path": "./tmp/repaired_ligand_data.csv",
        "index_cols": ["ligand_new_did"],
    },
    {
        "name": "l4",
        "csv_path": "./tmp/fragments.csv",
        "index_cols": ["fragment_DID"],
    },
    {
        "name": "l5",
        "csv_path": "./tmp/IRL_filtered.csv",
        "index_cols": ["DID"],
    },
]

# 可选调优：加载/查询期的 SQLite PRAGMA
PRAGMAS = [
    ("journal_mode", "OFF"),
    ("synchronous", "OFF"),
    ("temp_store", "MEMORY"),
    ("cache_size", "-200000"),   # 负数=KB；这里约200MB缓存
    ("mmap_size", str(1 << 30)), # 1GB mmap（系统允许的话）
]


# ======================
# 数据集封装（每个数据集一个独立的内存 SQLite）
# ======================
class Dataset:
    def __init__(self, name: str, csv_path: str, index_cols: List[str]):
        self.name = name
        self.csv_path = csv_path
        self.index_cols = index_cols
        self.headers: List[str] = []
        self.row_count: int = 0

        # 独立的内存库连接（常驻）
        self.conn: Optional[sqlite3.Connection] = None
        self._loaded: bool = False

    def _open_conn(self) -> sqlite3.Connection:
        conn = sqlite3.connect(":memory:", check_same_thread=False)
        conn.row_factory = sqlite3.Row
        for k, v in PRAGMAS:
            conn.execute(f"PRAGMA {k}={v}")
        return conn

    def load(self):
        if not os.path.exists(self.csv_path):
            raise FileNotFoundError(f"{self.csv_path} not found")

        t0 = time.time()
        self.conn = self._open_conn()
        cur = self.conn.cursor()

        # 读取表头，并创建全 TEXT 的表（最通用）
        with open(self.csv_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            self.headers = reader.fieldnames or []
            if not self.headers:
                raise RuntimeError(f"{self.csv_path} has no header")
            # 表名统一叫 t
            cols_sql = ", ".join([f'"{c}" TEXT' for c in self.headers])
            cur.execute(f'CREATE TABLE t ({cols_sql})')

            # 批量插入
            placeholders = ", ".join(["?"] * len(self.headers))
            insert_sql = f'INSERT INTO t VALUES ({placeholders})'

            batch = []
            batch_size = 1000
            n = 0
            for row in reader:
                n += 1
                batch.append([row.get(c, "") for c in self.headers])
                if len(batch) >= batch_size:
                    cur.executemany(insert_sql, batch)
                    batch.clear()
            if batch:
                cur.executemany(insert_sql, batch)

            self.row_count = n

        # 为索引列建索引（B-tree）
        for col in self.index_cols:
            if col not in self.headers:
                raise RuntimeError(f"Index column '{col}' not found in {self.csv_path}")
            cur.execute(f'CREATE INDEX idx_{col} ON t("{col}")')

        self.conn.commit()
        self._loaded = True
        t1 = time.time()
        print(f"[{self.name}] Loaded {self.row_count} rows from {self.csv_path} in {t1 - t0:.2f}s "
              f"(indexes: {', '.join(self.index_cols)})")

    def query_by_index(self, col: str, value: str, limit: Optional[int] = None) -> List[Dict[str, Any]]:
        if not self._loaded or self.conn is None:
            raise RuntimeError("dataset not loaded")

        if col not in self.headers:
            raise KeyError(f"column '{col}' not found in dataset '{self.name}'")

        q = f'SELECT * FROM t WHERE "{col}" = ?'
        params: List[Any] = [value]
        if limit is not None and limit > 0:
            q += " LIMIT ?"
            params.append(int(limit))

        rows = self.conn.execute(q, params).fetchall()
        return [dict(r) for r in rows]

    def stats(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "rows": self.row_count,
            "headers": self.headers,
            "index_cols": self.index_cols,
        }


# ======================
# 应用：加载 5 个数据集，暴露 API
# ======================
app = FastAPI(title="CSV in-memory query with SQLite (indexes on selected columns)")

# Allow CORS for any origin
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],        # Allow all origins
    allow_credentials=True,     # Allow cookies/Authorization headers
    allow_methods=["*"],        # Allow all HTTP methods
    allow_headers=["*"],        # Allow all headers
)

DATASETS: Dict[str, Dataset] = {
    cfg["name"]: Dataset(cfg["name"], cfg["csv_path"], cfg["index_cols"])
    for cfg in DATASETS_CONFIG
}

# Global variable to store combined stats in memory
COMBINED_STATS: Dict[str, Any] = {}


@app.on_event("startup")
def _startup():
    # 进程启动时一次性加载 5 个 CSV 到内存
    for name, ds in DATASETS.items():
        ds.load()
    
    # Load all stats files into memory
    _load_stats_files()


def _load_stats_files():
    """Load all statistics files into memory during startup"""
    global COMBINED_STATS
    
    stats_files = [
        "tmp/processing_stats.json",
        "tmp/repair_stats.json", 
        "tmp/step6_7_stats.json",
        "tmp/step8_stats.json",
        "tmp/step9_stats.json",
    ]
    
    COMBINED_STATS = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "available_stats": [],
        "processing_stats": None,
        "repair_stats": None,
        "step6_7_stats": None,
        "step8_stats": None,
        "step9_stats": None
    }
    
    for file_path in stats_files:
        try:
            if os.path.exists(file_path):
                with open(file_path, 'r', encoding='utf-8') as f:
                    data = ujson.load(f)
                    COMBINED_STATS["available_stats"].append(file_path)
                    
                    if "processing_stats.json" in file_path:
                        COMBINED_STATS["processing_stats"] = data
                    elif "repair_stats.json" in file_path:
                        COMBINED_STATS["repair_stats"] = data
                    elif "step6_7_stats.json" in file_path:
                        COMBINED_STATS["step6_7_stats"] = data
                    elif "step8_stats.json" in file_path:
                        COMBINED_STATS["step8_stats"] = data
                    elif "step9_stats.json" in file_path:
                        COMBINED_STATS["step9_stats"] = data
        except Exception as e:
            print(f"Warning: Could not load {file_path}: {e}")
    
    print(f"Loaded {len(COMBINED_STATS['available_stats'])} stats files into memory")


# ---- Pydantic models（便于文档） ----
class QueryResponse(BaseModel):
    dataset: str
    index_col: str
    value: str
    count: int
    rows: List[Dict[str, Any]]


# ---- API ----

@app.get("/datasets")
def list_datasets():
    return JSONResponse({name: ds.stats() for name, ds in DATASETS.items()})


@app.get("/query/{dataset}/{index_col}", response_model=QueryResponse)
def query(dataset: str,
          index_col: str,
          value: str = Query(..., description="等值匹配的索引值"),
          limit: Optional[int] = Query(None, gt=0, description="可选，限制返回条数")):
    ds = DATASETS.get(dataset)
    if not ds:
        raise HTTPException(status_code=404, detail=f"dataset '{dataset}' not found")
    try:
        rows = ds.query_by_index(index_col, value, limit=limit)
    except KeyError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {
        "dataset": dataset,
        "index_col": index_col,
        "value": value,
        "count": len(rows),
        "rows": rows,
    }


@app.get("/healthz")
def healthz():
    return {"ok": True}


@app.get("/api/db/stats")
def get_db_stats():
    """Return combined statistics from memory (loaded during startup)"""
    return JSONResponse(COMBINED_STATS)


# 本地启动：uvicorn app:app --reload
if __name__ == "__main__":
    import uvicorn
    # uvicorn.run("server:app", host="0.0.0.0", port=3041, reload=True)
    uvicorn.run("server:app", host="0.0.0.0", port=3041)
