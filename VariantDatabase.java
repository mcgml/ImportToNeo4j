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

    private Label patientLabel = DynamicLabel.label("Patient");
    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label variantLabel = DynamicLabel.label("Variant");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label geneLabel = DynamicLabel.label("Gene");
    private Label featureLabel = DynamicLabel.label("Feature");
    private Label codingSequenceVariantLabel = DynamicLabel.label("CodingSequenceVariant")
    private Label downstreamgeneVariantLabel = DynamicLabel.label(DownstreamgeneVariant)
    private Label featureElongationLabel = DynamicLabel.label(FeatureElongation)
    private Label featureTruncationLabel = DynamicLabel.label(FeatureTruncation)
    private Label fivePrimeUtrVariantLabel = DynamicLabel.label(FivePrimeUtrVariant)
    private Label frameshiftVariantLabel = DynamicLabel.label(FrameshiftVariant)
    private Label incompleteTerminalCodonVariantLabel = DynamicLabel.label(IncompleteTerminalCodonVariant)
    private Label inframeDeletionLabel = DynamicLabel.label(InframeDeletion)
    private Label inframeInsertionLabel = DynamicLabel.label(InframeInsertion)
    private Label initiatorCodonVariantLabel = DynamicLabel.label(InitiatorCodonVariant)
    private Label intergenicVariantLabel = DynamicLabel.label(IntergenicVariant)
    private Label intronVariantLabel = DynamicLabel.label(IntronVariant)
    private Label matureMiRNAVariantLabel = DynamicLabel.label(MatureMiRNAVariant)
    private Label missenseVariantLabel = DynamicLabel.label(MissenseVariant)
    private Label nmdTranscriptVariantLabel = DynamicLabel.label(NmdTranscriptVariant)
    private Label nonCodingTranscriptExonVariantLabel = DynamicLabel.label(NonCodingTranscriptExonVariant)
    private Label nonCodingTranscriptVariantLabel = DynamicLabel.label(NonCodingTranscriptVariant)
    private Label regulatoryRegionAblationLabel = DynamicLabel.label(RegulatoryRegionAblation)
    private Label regulatoryRegionAmplificationLabel = DynamicLabel.label(RegulatoryRegionAmplification )
    private Label regulatoryRegionVariantLabel = DynamicLabel.label(RegulatoryRegionVariant)
    private Label spliceAcceptorVariantLabel = DynamicLabel.label(SpliceAcceptorVariant)
    private Label spliceDonorVariantLabel = DynamicLabel.label(SpliceDonorVariant)
    private Label spliceRegionVariantLabel = DynamicLabel.label(SpliceRegionVariant)
    private Label stopGainedLabel = DynamicLabel.label(StopGained)
    private Label stopLostLabel = DynamicLabel.label(StopLost)
    private Label stopRetainedVariantLabel = DynamicLabel.label(StopRetainedVariant)
    private Label synonymousVariantLabel = DynamicLabel.label(SynonymousVariant)
    private Label tfBindingSiteVariantLabel = DynamicLabel.label(TfBindingSiteVariant)
    private Label tfbsAblationLabel = DynamicLabel.label(TfbsAblation)
    private Label tfbsAmplificationLabel = DynamicLabel.label(TfbsAmplification)
    private Label threePrimeUtrVariantLabel = DynamicLabel.label(ThreePrimeUtrVariant)
    private Label transcriptAblationLabel = DynamicLabel.label(TranscriptAblation)
    private Label transcriptAmplificationLabel = DynamicLabel.label(TranscriptAmplification)
    private Label upstreamGeneVariantLabel = DynamicLabel.label(UpstreamGeneVariant)






    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private ArrayList<VariantContext> variants = new ArrayList<>();
    private HashMap<String, Long> patientNodeIds = new HashMap<>(); //PatientID:NodeID
    private HashMap<String, Long> sampleNodeIds = new HashMap<>(); //sampleID:NodeID
    private HashMap<String, Long> variantNodeIds = new HashMap<>(); //VariantID:NodeID
    private HashMap<String, Long> annotationNodeIds = new HashMap<>(); //AnnotationID:NodeID
    private HashMap<String, Long> geneNodeIds = new HashMap<>(); //GeneID:NodeID
    private HashMap<String, Long> featureNodeIds = new HashMap<>(); //Feature:NodeID

    public VariantDatabase(VCFFileReader vcfFileReader, File neo4jDBPath){
        this.vcfFileReader = vcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_SAMPLE,
        HAS_VARIANT,
        IN_GENE,
        HAS_FEATURE,
        HAS_ANNOTATION
    }

    public void createDatabase(){

        log.log(Level.INFO, "Starting database ...");

        //create DB
        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }

    public void createIndexes() {
        Neo4j.createConstraint(graphDb, sampleLabel, "SampleID");
        Neo4j.createConstraint(graphDb, variantLabel, "VariantID");
        Neo4j.createConstraint(graphDb, annotationLabel, "AnnotationID");
        Neo4j.createConstraint(graphDb, geneLabel, "Symbol");
        Neo4j.createConstraint(graphDb, featureLabel, "Feature");
        Neo4j.createConstraint(graphDb, patientLabel, "PatientID");
    }

    public void addPatientNodes(){
        log.log(Level.INFO, "Adding patients nodes ...");
        //TODO
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

            }

        }

    }

    public void loadVariantsIntoMemory(){
        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (!variant.isFiltered() && !variant.isVariant()){
                variants.add(variant);
            }

        }
    }

    public void addVariants(){

        log.log(Level.INFO, "Adding variants ...");

        Node node;

        //read VCF file
        for (VariantContext variant : variants){

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

                    //add transcript node
                    if (!Neo4j.hasNode(graphDb, featureLabel, "Feature", annotation.getFeature())){

                        try (Transaction tx = graphDb.beginTx()) {

                            node = graphDb.createNode();
                            node.addLabel(featureLabel);

                            node.setProperty("Feature", annotation.getFeature());
                            node.setProperty("FeatureType", annotation.getFeatureType());

                            featureNodeIds.put(annotation.getFeature(), node.getId());

                            tx.success();
                        }

                    } else {
                        featureNodeIds.put(annotation.getFeature(), Neo4j.getUniqueMergedNode(graphDb, featureLabel, "Symbol", annotation.getSymbol()).getId());
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
                            if(annotation.getStrand() != null) node.setProperty("Strand", annotation.getStrand());

                            annotationNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature(), node.getId());

                            tx.success();
                        }

                        //link variant to annotation
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
