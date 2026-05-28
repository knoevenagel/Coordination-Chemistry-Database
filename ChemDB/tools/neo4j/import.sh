export NEO4J_IMAGE=neo4j:2025.07
export DB_NAME=neo4j

docker compose down

docker run --rm \
  -v "../../tmp":/import \
  -v "$PWD/neo4j/data":/data \
  -v "$PWD/neo4j/logs":/logs \
  $NEO4J_IMAGE \
  neo4j-admin database import full "$DB_NAME" \
    --verbose \
    --overwrite-destination=true \
    --nodes=/import/neo4j_metals.csv \
    --nodes=/import/neo4j_complexes.csv \
    --nodes=/import/neo4j_ligands.csv \
    --nodes=/import/neo4j_repaired_ligands.csv \
    --nodes=/import/neo4j_fragments.csv \
    --nodes=/import/neo4j_irl.csv \
    --relationships=/import/neo4j_m_l1_relationships.csv \
    --relationships=/import/neo4j_l1_l3_relationships.csv \
    --relationships=/import/neo4j_l2_l3_relationships.csv \
    --relationships=/import/neo4j_l3_l4_relationships.csv \
    --relationships=/import/neo4j_l4_l5_relationships.csv \
    --delimiter="," \
    --array-delimiter="|" \
    --trim-strings=true \
    --bad-tolerance=1000
