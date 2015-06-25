package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.File;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 23/06/15.
 */
public class VariantDatabase {

    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label variantLabel = DynamicLabel.label("Variant");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label geneLabel = DynamicLabel.label("Gene");
    private Label phenotypeLabel = DynamicLabel.label("Phenotype");

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private HashMap<String, Long> sampleNodeIds = new HashMap<>(); //sampleID:NodeID
    private HashMap<String, Long> variantNodeIds = new HashMap<>(); //VariantID:NodeID
    private HashMap<String, Long> annotationNodeIds = new HashMap<>(); //AnnotationID:NodeID
    private HashMap<String, Long> geneNodeIds = new HashMap<>(); //GeneID:NodeID
    private HashMap<String, Long> phenotypeNodeIds = new HashMap<>(); //PhenotypeID:NodeID

    //TODO add allele frequencies
    //TODO add custom variant classification

    public VariantDatabase(VCFFileReader vcfFileReader, File neo4jDBPath){
        this.vcfFileReader = vcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_VARIANT,
        IN_GENE,
        HAS_ANNOTATION,
        HAS_PHENOTYPE,
        HAS_SAMPLE,
        ASSOCIATED_PHENOTYPE
    }

    public void createDatabase(){

        log.log(Level.INFO, "Starting database ...");

        //create DB
        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }

    public void createIndexes(){
        Neo4j.createConstraint(graphDb, sampleLabel, "SampleID");
        Neo4j.createConstraint(graphDb, variantLabel, "VariantID");
        Neo4j.createConstraint(graphDb, annotationLabel, "AnnotationID");
        Neo4j.createConstraint(graphDb, geneLabel, "Symbol");
    }

    public void addSampleNodes(){

        log.log(Level.INFO, "Adding sample nodes ...");

        Node node;

        for (String sampleID : vcfFileReader.getFileHeader().getSampleNamesInOrder()){

            if (!Neo4j.hasNode(graphDb, sampleLabel, "SampleID", sampleID)){

                //add sampleID
                try (Transaction tx = graphDb.beginTx()) {

                    node = graphDb.createNode();
                    node.addLabel(sampleLabel);
                    node.setProperty("SampleID", sampleID);

                    //store node Id for later
                    sampleNodeIds.put(sampleID, node.getId());

                    tx.success();
                }

            } else {
                log.log(Level.WARNING, sampleID + " already exists in database.");
                //TODO: delete sample and relationships?
            }

        }

    }

    public void addVariants(){

        log.log(Level.INFO, "Adding variants ...");

        Node node;
        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //loop over alleles
            for (Allele allele : variant.getAlleles()){

                //skip wildtype alleles
                if (allele.isReference()) continue;

                if (!Neo4j.hasNode(graphDb, variantLabel, "VariantID",  variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString())){

                    //add variant
                    try (Transaction tx = graphDb.beginTx()) {

                        node = graphDb.createNode();
                        node.addLabel(variantLabel);

                        node.setProperty("VariantID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString());
                        node.setProperty("Contig", variant.getContig());
                        node.setProperty("Position", variant.getStart());
                        node.setProperty("Reference", variant.getReference().getBaseString());
                        node.setProperty("Alternative", allele.getBaseString());

                        //store node Id for later
                        variantNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString(), node.getId());

                        tx.success();
                    }

                } else {
                    //store node Id for later
                    variantNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString(),
                            Neo4j.getUniqueMergedNode(graphDb, variantLabel, "VariantID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString()).getId());
                }

            }

        }

    }

    public void addFunctionalAnnotations(){

        log.log(Level.INFO, "Adding functional annotations ...");

        Node node;
        HashMap<String, HashSet<VEPAnnotation>> splitAnnotations = new HashMap<>();
        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotation splitAnnotation = new VEPAnnotation((String)variant.getAttribute("CSQ"));
                splitAnnotation.parseAnnotation();

                if (!splitAnnotations.containsKey(splitAnnotation.getAllele())) splitAnnotations.put(splitAnnotation.getAllele(), new HashSet<VEPAnnotation>());
                splitAnnotations.get(splitAnnotation.getAllele()).add(splitAnnotation);

            } catch (ClassCastException e){

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")){

                    //multiple annotations
                    VEPAnnotation splitAnnotation = new VEPAnnotation(annotation);
                    splitAnnotation.parseAnnotation();

                    if (!splitAnnotations.containsKey(splitAnnotation.getAllele())) splitAnnotations.put(splitAnnotation.getAllele(), new HashSet<VEPAnnotation>());
                    splitAnnotations.get(splitAnnotation.getAllele()).add(splitAnnotation);

                }

            }

            //loop over alleles
            for (Allele allele : variant.getAlleles()){

                //skip wildtype alleles
                if (allele.isReference()) continue;

                //skip missing alleles
                if (!splitAnnotations.containsKey(allele.getBaseString())) {
                    log.log(Level.WARNING, "Missing annotation for: " + variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString());
                    continue;
                }

                //loop over annotations by for this allele
                for (VEPAnnotation annotation : splitAnnotations.get(allele.getBaseString())){

                    //transcript, ensg, symbol
                    if (annotation.getFeature() == null || annotation.getGene() == null || annotation.getSymbol() == null) continue;

                    //add gene node
                    if (!Neo4j.hasNode(graphDb, geneLabel, "Symbol", annotation.getSymbol())){

                        try (Transaction tx = graphDb.beginTx()) {

                            node = graphDb.createNode();
                            node.addLabel(geneLabel);

                            node.setProperty("GeneID", annotation.getGene());
                            node.setProperty("Symbol", annotation.getSymbol());

                            geneNodeIds.put(annotation.getGene(), node.getId());

                            tx.success();
                        }

                    } else {
                        geneNodeIds.put(annotation.getSymbol(), Neo4j.getUniqueMergedNode(graphDb, geneLabel, "Symbol", annotation.getSymbol()).getId());
                    }

                    //make annotation node
                    if (!Neo4j.hasNode(graphDb, annotationLabel, "AnnotationID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature())){

                        try (Transaction tx = graphDb.beginTx()) {

                            node = graphDb.createNode();
                            node.addLabel(annotationLabel);
                            node.setProperty("AnnotationID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature());

                            if(annotation.getExon() != null) {
                                String[] fields = annotation.getExon().split("/");
                                node.setProperty("Exon", fields[0]);
                                node.setProperty("TotalExons", fields[1]);
                            }
                            if(annotation.getFeature() != null) node.setProperty("Feature", annotation.getFeature());
                            if(annotation.getFeatureType() != null) node.setProperty("FeatureType", annotation.getFeatureType());
                            if(annotation.getHgvsCoding() != null) node.setProperty("HGVSc", annotation.getHgvsCoding());
                            if(annotation.getHgvsProtein() != null) node.setProperty("HGVSp", annotation.getHgvsProtein());
                            if(annotation.getIntron() != null) {
                                String[] fields = annotation.getIntron().split("/");
                                node.setProperty("Intron", fields[0]);
                                node.setProperty("TotalIntrons", fields[1]);
                            }
                            if(annotation.getPolyphen() != null) node.setProperty("Polyphen", annotation.getPolyphen());
                            if(annotation.getSift() != null) node.setProperty("Sift", annotation.getSift());
                            if(annotation.getStrand() != null) node.setProperty("Strand", annotation.getStrand());


                            annotationNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature(), node.getId());

                            tx.success();
                        }

                        //link annotation to variant
                        //todo store consequence
                        try (Transaction tx = graphDb.beginTx()) {

                            Relationship relationship = graphDb.getNodeById(variantNodeIds.get(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString()))
                                    .createRelationshipTo(node, relTypes.HAS_ANNOTATION);

                            tx.success();

                        } //done creating relationship

                    } else {
                        annotationNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature(),
                                Neo4j.getUniqueMergedNode(graphDb, annotationLabel, "AnnotationID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature()).getId());
                    }

                } //done looping over annotations

            } //done looping over alleles

            splitAnnotations.clear();

        } //done reading VCF

    }

    //TODO check why? duplicate relationships
    public void linkSamplesToVariants(){

        log.log(Level.INFO, "Linking samples to variants ...");

        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()){
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                //loop over alleles
                for (Allele allele : genotype.getAlleles()){

                    //skip wildtype alleles
                    if (allele.isReference()) continue;

                    //create relationship
                    try (Transaction tx = graphDb.beginTx()) {

                        Relationship relationship = graphDb.getNodeById(sampleNodeIds.get(genotype.getSampleName()))
                                .createRelationshipTo(graphDb.getNodeById(variantNodeIds.get(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString())), relTypes.HAS_VARIANT);

                        if (genotype.isHet()) relationship.setProperty("Het", true); else relationship.setProperty("Het", false);
                        if (genotype.isHetNonRef()) relationship.setProperty("HetNonRef", true); else relationship.setProperty("HetNonRef", false);
                        if (genotype.isHomVar()) relationship.setProperty("HomVar", true); else relationship.setProperty("HomVar", false);
                        if (genotype.isMixed()) relationship.setProperty("Mixed", true); else relationship.setProperty("Mixed", false);
                        relationship.setProperty("GQ", genotype.getGQ());

                        tx.success();

                    }

                }
            }

        }

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");

        graphDb.shutdown();
    }

}
