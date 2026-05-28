# proxy.py
import typing as T
from fastapi import FastAPI, Request, Response
from fastapi.responses import StreamingResponse, JSONResponse
import httpx
from starlette.background import BackgroundTask


app = FastAPI(title="Simple Reverse Proxy")

BACKENDS = {
    "SIM": "http://localhost:3042",
    "QRY": "http://localhost:3041",
    "SPA": "http://localhost:3043",
    "EMB": "http://localhost:3044",
}

HOP_BY_HOP = {
    "connection", "keep-alive", "proxy-authenticate", "proxy-authorization",
    "te", "trailers", "transfer-encoding", "upgrade",
}

client = httpx.AsyncClient(follow_redirects=False, timeout=None)

def choose_backend(path: str) -> str:
    # 优先级顺序很重要：更具体的规则要先匹配
    if path.startswith("/api/embeddings"):
        return BACKENDS["EMB"]
    if path.startswith("/api/similarity"):
        return BACKENDS["SIM"]
    if path.startswith("/query/"):
        return BACKENDS["QRY"]
    if path.startswith("/api/db/"):
        return BACKENDS["QRY"]
    return BACKENDS["SPA"]

def build_target_url(base: str, req: Request) -> str:
    # 路径保持不变，查询串原样透传
    if req.url.query:
        return f"{base}{req.url.path}?{req.url.query}"
    return f"{base}{req.url.path}"

def filter_request_headers(headers: T.Mapping[str, str]) -> dict:
    new = {}
    for k, v in headers.items():
        lk = k.lower()
        if lk in HOP_BY_HOP:
            continue
        if lk == "host":
            continue
        # 接收端自己决定压缩，避免端到端编码问题
        if lk == "accept-encoding":
            continue
        new[k] = v
    return new

def filter_response_headers(headers: T.Mapping[str, str]) -> dict:
    new = {}
    for k, v in headers.items():
        if k.lower() in HOP_BY_HOP:
            continue
        # 让 ASGI 框架处理压缩/长度
        if k.lower() in {"content-encoding", "content-length"}:
            continue
        new[k] = v
    return new

async def proxy(req: Request) -> Response:
    base = choose_backend(req.url.path)
    url = build_target_url(base, req)
    method = req.method.upper()
    headers = filter_request_headers(req.headers)

    # 取请求体（支持任意方法）
    body = await req.body()

    req_up = client.build_request(method, url, headers=headers, content=body)
    r = await client.send(req_up, stream=True)

    resp_headers = filter_response_headers(r.headers)

    async def iter_bytes():
        async for chunk in r.aiter_bytes():
            yield chunk

    # 确保流结束后正确关闭上游连接
    return StreamingResponse(
        iter_bytes(),
        status_code=r.status_code,
        headers=resp_headers,
        background=BackgroundTask(r.aclose),
    )


# 明确注册所有常见方法，并在末尾放一个兜底的 catch-all
for _method in ["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS", "HEAD"]:
    app.add_api_route("/{full_path:path}", proxy, methods=[_method])

@app.get("/healthz")
async def healthz():
    return {"ok": True}

@app.on_event("shutdown")
async def _shutdown():
    await client.aclose()


# 本地启动：uvicorn app:app --reload
if __name__ == "__main__":
    import uvicorn
    # uvicorn.run("proxy:app", host="0.0.0.0", port=3040, reload=True)
    uvicorn.run("proxy:app", host="0.0.0.0", port=3040)
