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
    public static Node addNode(final GraphDatabaseService graphDb, final Label label, HashMap<String, Object> properties){

        Node node;

        try ( Transaction tx = graphDb.beginTx() )
        {
            node = graphDb.createNode( label );

            for (Map.Entry<String, Object> property : properties.entrySet()){
                node.setProperty(property.getKey(), property.getValue());
            }

            tx.success();
        }

        return node;

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
    public static Node matchOrCreateUniqueNode(final GraphDatabaseService graphDb, Label label, String field, Object value) throws InvalidPropertiesFormatException{
        ArrayList<Node> nodes = getNodes(graphDb, label, field, value);

        if (nodes.size() == 1){
            return nodes.get(0);
        } else if (nodes.size() > 1){
            throw new InvalidPropertiesFormatException("Multiple nodes already exists for: " + label.name() + " " + field + " " + value);
        }

        HashMap<String, Object> props = new HashMap<>();
        props.put(field, value);

        return addNode(graphDb, label, props);
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
    public static void createRelationship(final GraphDatabaseService graphDb, Node node1, Node node2, RelationshipType type, HashMap<String, Object> properties, boolean checkDuplicates){

        if (!checkDuplicates){
            if (hasRelationship(graphDb, node1, node2, type, Direction.OUTGOING)){
                return;
            }
        }

        //add relationship
        try (Transaction tx = graphDb.beginTx()) {

            Relationship relationship = node1.createRelationshipTo(node2, type);

            //set properties
            for (Map.Entry<String, Object> property : properties.entrySet()){
                relationship.setProperty(property.getKey(), property.getValue());
            }

            tx.success();
        }

    }
    public static boolean hasRelationship(final GraphDatabaseService graphDb, Node node1, Node node2, RelationshipType type, Direction direction){

        //check if relationship already exists
        try ( Transaction tx = graphDb.beginTx() ){

            for (Relationship relationship : node1.getRelationships(type, direction)){

                if (relationship.getOtherNode(node1).getId() == node2.getId()){
                    return true;
                }

            }

        }

        return false;
    }
    public static void addNodeProperty(final GraphDatabaseService graphDb, Node node, HashMap<String, Object> properties){

        try (Transaction tx = graphDb.beginTx()) {

            //set properties
            for (Map.Entry<String, Object> property : properties.entrySet()){
                node.setProperty(property.getKey(), property.getValue());
            }

            tx.success();
        }

    }
    public static void addNodeLabel(final GraphDatabaseService graphDb, Node node, Label label){

        try (Transaction tx = graphDb.beginTx()) {

            //add label
            node.addLabel(label);

            tx.success();
        }

    }
    public static ArrayList<Node> findNeighbourNodes(final GraphDatabaseService graphDb, Node startNode, Label endLabel, Direction direction){

        Node tempNode;
        ArrayList<Node> nodes = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx() ){

            for (Relationship relationship : startNode.getRelationships(direction)){

                //get connecting node
                tempNode = relationship.getOtherNode(startNode);

                //check node has required label
                if (tempNode.hasLabel(endLabel)){
                    nodes.add(tempNode);
                }

            }

        }

        return nodes;
    }
    public static boolean isNeighbourNodeWithSuppliedProperties(final GraphDatabaseService graphDb, Node startNode, Node endNode, Direction direction, RelationshipType relationshipType, HashMap<String, Object> properties){

        boolean allPropertiesMatched;

        try ( Transaction tx = graphDb.beginTx() ){
            for (Relationship relationship : startNode.getRelationships(direction, relationshipType)){

                //check if node is endnode
                if (relationship.getOtherNode(startNode).equals(endNode)){
                    allPropertiesMatched = true;

                    for (Map.Entry<String, Object> property : properties.entrySet()){

                        if (!relationship.getProperty(property.getKey()).equals(property.getValue())){
                            allPropertiesMatched = false;
                            break;
                        }

                    }

                    if (allPropertiesMatched){
                        return true;
                    }
                }

            }
        }

        return false;
    }
}