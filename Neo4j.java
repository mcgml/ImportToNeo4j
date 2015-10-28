package nhs.genetics.cardiff;

import org.neo4j.graphdb.*;
import org.neo4j.graphdb.schema.IndexDefinition;
import org.neo4j.graphdb.schema.Schema;

import java.util.*;
import java.util.concurrent.TimeUnit;

public class Neo4j{

    //labels
    private static Label sampleLabel = DynamicLabel.label("Sample");
    private static Label variantLabel = DynamicLabel.label("Variant");
    private static Label autoChromosomeLabel = DynamicLabel.label("Autosome");
    private static Label xChromosomeLabel = DynamicLabel.label("X");
    private static Label yChromosomeLabel = DynamicLabel.label("Y");
    private static Label mtChromosomeLabel = DynamicLabel.label("MT");
    private static Label annotationLabel = DynamicLabel.label("Annotation");
    private static Label symbolLabel = DynamicLabel.label("Symbol");
    private static Label canonicalLabel = DynamicLabel.label("Canonical");
    private static Label featureLabel = DynamicLabel.label("Feature");
    private static Label runInfoLabel = DynamicLabel.label("RunInfo");
    private static Label virtualPanelLabel = DynamicLabel.label("VirtualPanel");
    private static Label userLabel = DynamicLabel.label("User");

    //relationships
    private static RelationshipType hasHetVariantRelationship = DynamicRelationshipType.withName("HAS_HET_VARIANT");
    private static RelationshipType hasHomVariantRelationship = DynamicRelationshipType.withName("HAS_HOM_VARIANT");
    private static RelationshipType hasInSymbolRelationship = DynamicRelationshipType.withName("IN_SYMBOL");
    private static RelationshipType hasInFeatureRelationship = DynamicRelationshipType.withName("IN_FEATURE");
    private static RelationshipType hasUnknownConsequenceRelationship = DynamicRelationshipType.withName("HAS_UNKNOWN_CONSEQUENCE");
    private static RelationshipType hasAnalysisRelationship = DynamicRelationshipType.withName("HAS_ANALYSIS");
    private static RelationshipType hasDesignedByRelationship = DynamicRelationshipType.withName("DESIGNED_BY");
    private static RelationshipType hasContainsSymbolRelationship = DynamicRelationshipType.withName("CONTAINS_SYMBOL");
    private static RelationshipType hasProteinCodingBiotypeRelationship = DynamicRelationshipType.withName("HAS_PROTEIN_CODING_BIOTYPE");
    private static RelationshipType hasAssignedPathogenicityRelationship = DynamicRelationshipType.withName("HAS_ASSIGNED_PATHOGENICITY");
    private static RelationshipType hasUserCommentRelationship = DynamicRelationshipType.withName("HAS_USER_COMMENT");

    //population frequencies
    public enum Population {
        onekGPhase3_EAS_AF,
        onekGPhase3_EUR_AF,
        onekGPhase3_AFR_AF,
        onekGPhase3_AMR_AF,
        onekGPhase3_SAS_AF,
        ExAC_AFR_AF,
        ExAC_AMR_AF,
        ExAC_EAS_AF,
        ExAC_FIN_AF,
        ExAC_NFE_AF,
        ExAC_OTH_AF,
        ExAC_SAS_AF
    }

    public static void shutdownDatabase(final GraphDatabaseService graphDb){
        graphDb.shutdown();
    }
    public static void registerShutdownHook( final GraphDatabaseService graphDb )
    {
        // Registers a shutdown hook for the Neo4j instance so that it
        // shuts down nicely when the VM exits (even if you "Ctrl-C" the
        // running application).
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
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
            tx.success();
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

            tx.success();
        }

        return nodes;
    }
    public static ArrayList<Long> getNodeIds(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        Node node;
        ArrayList<Long> nodeIds = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx() )
        {
            try ( ResourceIterator<Node> users = graphDb.findNodes( label, field, value ) )
            {

                while ( users.hasNext() )
                {
                    node = users.next();
                    nodeIds.add(node.getId());
                }

            }

            tx.success();
        }

        return nodeIds;
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
              Result result = graphDb.execute(cypherQuery) )
        {
            while ( result.hasNext() )
            {
                results.add(result.next());
            }

            tx.success();
        }

        return results;
    }
    public static void createRelationship(final GraphDatabaseService graphDb, Node node1, Node node2, RelationshipType type, HashMap<String, Object> properties){

        if (hasRelationship(graphDb, node1, node2, type, Direction.OUTGOING)){
            return;
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

            tx.success();
        }

        return false;
    }
    public static void addNodeProperties(final GraphDatabaseService graphDb, Node node, HashMap<String, Object> properties){

        try (Transaction tx = graphDb.beginTx()) {

            //set properties
            for (Map.Entry<String, Object> property : properties.entrySet()){
                if (!node.hasProperty(property.getKey())) node.setProperty(property.getKey(), property.getValue());
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

            tx.success();
        }

        return nodes;
    }
    public static ArrayList<Node> findNeighbourNodesWithParameters(final GraphDatabaseService graphDb, Node startNode, Label endLabel, Direction direction, RelationshipType relationshipType, HashMap<String, Object> properties){

        ArrayList<Node> nodes = new ArrayList<>();

        try (Transaction tx = graphDb.beginTx()){
            for (Relationship relationship : startNode.getRelationships(direction, relationshipType)){
                boolean hasSameProperties = true;

                //get connecting node
                Node tempNode = relationship.getOtherNode(startNode);

                //check node has required label
                if (tempNode.hasLabel(endLabel)){

                    //check end node has properties
                    for (Map.Entry<String, Object> property : properties.entrySet()){

                        if (!tempNode.getProperty(property.getKey()).equals(property.getValue())){
                            hasSameProperties = false;
                            break;
                        }
                    }

                    if (hasSameProperties){
                        nodes.add(tempNode);
                    }

                }

            }
            tx.success();
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
            tx.success();
        }

        return false;
    }

    public static Label getSampleLabel() {
        return sampleLabel;
    }
    public static Label getVariantLabel() {
        return variantLabel;
    }
    public static Label getAutoChromosomeLabel() {
        return autoChromosomeLabel;
    }
    public static Label getXChromosomeLabel() {
        return xChromosomeLabel;
    }
    public static Label getYChromosomeLabel() {
        return yChromosomeLabel;
    }
    public static Label getMtChromosomeLabel() {
        return mtChromosomeLabel;
    }
    public static Label getAnnotationLabel() {
        return annotationLabel;
    }
    public static Label getSymbolLabel() {
        return symbolLabel;
    }
    public static Label getCanonicalLabel() {
        return canonicalLabel;
    }
    public static Label getFeatureLabel() {
        return featureLabel;
    }
    public static Label getRunInfoLabel() {
        return runInfoLabel;
    }
    public static Label getVirtualPanelLabel() {
        return virtualPanelLabel;
    }
    public static Label getUserLabel() {
        return userLabel;
    }
    public static RelationshipType getHasHetVariantRelationship() {
        return hasHetVariantRelationship;
    }
    public static RelationshipType getHasHomVariantRelationship() {
        return hasHomVariantRelationship;
    }
    public static RelationshipType getHasInSymbolRelationship() {
        return hasInSymbolRelationship;
    }
    public static RelationshipType getHasInFeatureRelationship() {
        return hasInFeatureRelationship;
    }
    public static RelationshipType getHasUnknownConsequenceRelationship() {
        return hasUnknownConsequenceRelationship;
    }
    public static RelationshipType getHasAnalysisRelationship() {
        return hasAnalysisRelationship;
    }
    public static RelationshipType getHasContainsSymbol() {
        return hasContainsSymbolRelationship;
    }
    public static RelationshipType getHasDesignedBy() {
        return hasDesignedByRelationship;
    }
    public static RelationshipType getHasProteinCodingBiotypeRel() {
        return hasProteinCodingBiotypeRelationship;
    }
    public static RelationshipType getHasAssignedPathogenicityRelationship() {
        return hasAssignedPathogenicityRelationship;
    }
    public static RelationshipType getHasUserCommentRelationship() {
        return hasUserCommentRelationship;
    }
    public static RelationshipType getHasProteinCodingBiotypeRelationship() {
        return hasProteinCodingBiotypeRelationship;
    }
    public static RelationshipType getHasContainsSymbolRelationship() {
        return hasContainsSymbolRelationship;
    }
    public static RelationshipType getHasDesignedByRelationship() {
        return hasDesignedByRelationship;
    }

}