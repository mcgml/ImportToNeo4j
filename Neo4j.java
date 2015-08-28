package nhs.genetics.cardiff;

import org.neo4j.graphdb.*;
import org.neo4j.graphdb.schema.IndexDefinition;
import org.neo4j.graphdb.schema.Schema;

import java.util.*;
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
    public static void addNode(final GraphDatabaseService graphDb, final Label label, HashMap<String, Object> properties){

        try ( Transaction tx = graphDb.beginTx() )
        {
            Node userNode = graphDb.createNode( label );

            for(Map.Entry<String, Object> property : properties.entrySet()){
                userNode.setProperty(property.getKey(), property.getValue());
            }

            tx.success();
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
    public static ArrayList<Map<String, Object>> runCypherQuery(final GraphDatabaseService graphDb, String cypherQuery){
        ArrayList<Map<String, Object>> results = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx();
              Result result = graphDb.execute( cypherQuery ) )
        {
            while ( result.hasNext() )
            {
                results.add(result.next());
            }

            tx.success();
        }

        return results;
    }
    public static void createRelationship(final GraphDatabaseService graphDb, Long nodeId1, Long nodeId2, RelationshipType type, HashMap<String, Object> properties){

        try (Transaction tx = graphDb.beginTx()) {

            Node node1 = graphDb.getNodeById(nodeId1);
            Node node2 = graphDb.getNodeById(nodeId2);

            Relationship relationship = node1.createRelationshipTo(node2, type);

            //set properties
            for (Map.Entry<String, Object> property : properties.entrySet()){
                relationship.setProperty(property.getKey(), property.getValue());
            }

            tx.success();
        }

    }
    public static void addNodeProperty(final GraphDatabaseService graphDb, Long nodeId, HashMap<String, Object> properties){

        try (Transaction tx = graphDb.beginTx()) {

            Node node = graphDb.getNodeById(nodeId);

            //set properties
            for (Map.Entry<String, Object> property : properties.entrySet()){
                node.setProperty(property.getKey(), property.getValue());
            }

            tx.success();
        }

    }
    public static void addNodeLabel(final GraphDatabaseService graphDb, Long nodeId, Label label){

        try (Transaction tx = graphDb.beginTx()) {

            Node node = graphDb.getNodeById(nodeId);

            //add label
            node.addLabel(label);

            tx.success();
        }

    }
    public static long matchOrCreateUniqueNode(final GraphDatabaseService graphDb, final Label label, final String field, final Object value) throws InvalidPropertiesFormatException{

        ArrayList<Node> nodes = new ArrayList<>();

        //check if node exists
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

        //add node
        if (nodes.size() == 0){
            try ( Transaction tx = graphDb.beginTx() )
            {
                Node userNode = graphDb.createNode( label );
                userNode.setProperty(field, value);
                nodes.add(userNode);

                tx.success();
            }
        }

        if (nodes.size() > 1){
            throw new InvalidPropertiesFormatException("Multiple nodes were present in DB for: " + label.name() + " " + field + " " + value);
        }

        return nodes.get(0).getId();
    }

}