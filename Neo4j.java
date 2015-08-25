package nhs.genetics.cardiff;

import org.neo4j.graphdb.*;
import org.neo4j.graphdb.schema.IndexDefinition;
import org.neo4j.graphdb.schema.Schema;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

public class Neo4j{
    public static void registerShutdownHook( final GraphDatabaseService graphDb )
    {
        // Registers a shutdown hook for the Neo4j instance so that it
        // shuts down nicely when the VM exits (even if you "Ctrl-C" the
        // running application).
        Runtime.getRuntime().addShutdownHook( new Thread()
        {
            @Override
            public void run()
            {
                graphDb.shutdown();
            }
        } );
    }
    public static void createIndexAndWait(final GraphDatabaseService graphDb, final Label label, final String property){

        IndexDefinition indexDefinition;

        try ( Transaction tx = graphDb.beginTx() )
        {
            Schema schema = graphDb.schema();
            indexDefinition = schema.indexFor( label )
                    .on(property)
                    .create();
            tx.success();
        }

        try ( Transaction tx = graphDb.beginTx() )
        {
            Schema schema = graphDb.schema();
            schema.awaitIndexOnline(indexDefinition, 10, TimeUnit.SECONDS);
        }

    }
    public static void createIndex(final GraphDatabaseService graphDb, final Label label, final String property){

        try ( Transaction tx = graphDb.beginTx() )
        {
            graphDb.schema()
                    .indexFor(label)
                    .on(property)
                    .create();

            tx.success();
        }

    }
    public static void createConstraint(final GraphDatabaseService graphDb, final Label label, final String property) {

        //TODO check constraint does not already exist

        try ( Transaction tx = graphDb.beginTx() )
        {
            graphDb.schema()
                    .constraintFor(label)
                    .assertPropertyIsUnique(property)
                    .create();

            tx.success();
        }

    }
    public static void dropIndex(final GraphDatabaseService graphDb, final Label label){

        try ( Transaction tx = graphDb.beginTx() )
        {
            for ( IndexDefinition indexDefinition : graphDb.schema()
                    .getIndexes(label) )
            {
                indexDefinition.drop();

            }

            tx.success();
        }

    }
    public static void addNode(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        try ( Transaction tx = graphDb.beginTx() )
        {
            Node userNode = graphDb.createNode( label );
            userNode.setProperty(field, value);

            tx.success();
        }

    }
    public static boolean hasNode(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        try ( Transaction tx = graphDb.beginTx() )
        {
            if (graphDb.findNode(label, field, value) == null){
                return false;
            } else {
                return true;
            }
        }

    }
    public static ArrayList<Node> getNodes(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        ArrayList<Node> nodes = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx() )
        {
            try ( ResourceIterator<Node> users = graphDb.findNodes( label, field, value ) )
            {

                while ( users.hasNext() )
                {
                    nodes.add( users.next() );
                }

            }

        }

        return nodes;
    }
    public static ArrayList<Long> getNodeIDs(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        Node node;
        ArrayList<Long> nodeIDs = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx() )
        {
            try ( ResourceIterator<Node> users = graphDb.findNodes( label, field, value ) )
            {

                while ( users.hasNext() )
                {
                    node = users.next();
                    nodeIDs.add(node.getId());
                }

            }

        }

        return nodeIDs;
    }
    public static Node getUniqueMergedNode(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        Node result = null;
        ResourceIterator<Node> resultIterator = null;
        try ( Transaction tx = graphDb.beginTx() )
        {
            HashMap<String, Object> parameters = new HashMap<>();

            String queryString = "MERGE (n:label {field: {value}}) RETURN n";
            parameters.put( "label", label.name() );
            parameters.put( "field", field );
            parameters.put( "value", value );

            resultIterator = graphDb.execute( queryString, parameters ).columnAs( "n" );
            result = resultIterator.next();

            tx.success();

            return result;
        }

    }
    public static ArrayList<Map<String, Object>> runCypherQuery(final GraphDatabaseService graphDb, String cypherQuery){

        ArrayList<Map<String, Object>> results = new ArrayList<>();

        try ( Transaction ignored = graphDb.beginTx();
              Result result = graphDb.execute( cypherQuery ) )
        {
            while ( result.hasNext() )
            {
                results.add(result.next());
            }
        }

        return results;

    }
    public static void createRelationship(final GraphDatabaseService graphDb, Long nodeId1, Long nodeId2, RelationshipType type){

        try (Transaction tx = graphDb.beginTx()) {

            Node node1 = graphDb.getNodeById(nodeId1);
            Node node2 = graphDb.getNodeById(nodeId2);

            Relationship relationship = node1.createRelationshipTo(node2, type);
            tx.success();
        }

    }
}