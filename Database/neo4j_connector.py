from neo4j import GraphDatabase

class Neo4jConnector:
    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    def test_connection(self):
        with self.driver.session() as session:
            result = session.run("RETURN 'Neo4j connection successful!' AS message")
            for record in result:
                print(record["message"])

# 修改为你的 Neo4j 配置
uri = "bolt://localhost:7687"  # 或 "neo4j://localhost:7687"（根据 Neo4j 版本而定）
username = "neo4j"
password = "zwl020115"

# 建立连接并测试
neo4j_conn = Neo4jConnector(uri, username, password)
neo4j_conn.test_connection()
neo4j_conn.close()