package nhs.genetics.cardiff;

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

    //TODO store everything as node not nodeId

    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader variantVcfFileReader, oneKgP3VcfFileReader;
    private ArrayList<VariantContext> variants = new ArrayList<>();
    private HashMap<String, HashSet<VEPAnnotation>> vepAnnotations = new HashMap<>(); //vepVar:ann
    private HashMap<GenomeVariant, Long> variantNodeIds = new HashMap<>();
    private HashMap<String, Long> patientNodeIds = new HashMap<>(); //PatientID:NodeID
    private HashMap<String, Long> sampleNodeIds = new HashMap<>(); //sampleID:NodeID
    private HashMap<GenomeVariant, HashMap<String, Long>> annotationNodeIds = new HashMap<>(); //genomeVariant:transcriptId=nodeId
    private HashMap<String, Long> symbolNodeIds = new HashMap<>(); //GeneID:NodeID
    private HashMap<String, Long> featureNodeIds = new HashMap<>(); //Feature:NodeID
    private String libraryId = "150716_D00501_0047_BHB092ADXX";

    private Label patientLabel = DynamicLabel.label("Patient");
    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label symbolLabel = DynamicLabel.label("Symbol");
    private Label canonicalLabel = DynamicLabel.label("Canonical");
    private Label featureLabel = DynamicLabel.label("Feature");
    private Label snpLabel = DynamicLabel.label("Snp");
    private Label indelLabel = DynamicLabel.label("Indel");

    public VariantDatabase(VCFFileReader variantVcfFileReader, VCFFileReader oneKgP3VcfFileReader, File neo4jDBPath){
        this.variantVcfFileReader = variantVcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
        this.oneKgP3VcfFileReader = oneKgP3VcfFileReader;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_SAMPLE,
        HAS_HET_VARIANT,
        HAS_HOM_VARIANT,
        IN_SYMBOL,
        IN_FEATURE,
        HAS_TRANSCRIPT_ABLATION_CONSEQUENCE,
        HAS_SPLICE_ACCEPTOR_VARIANT_CONSEQUENCE,
        HAS_SPLICE_DONOR_VARIANT_CONSEQUENCE,
        HAS_STOP_GAINED_CONSEQUENCE,
        HAS_FRAMESHIFT_VARIANT_CONSEQUENCE,
        HAS_STOP_LOST_CONSEQUENCE,
        HAS_START_LOST_CONSEQUENCE,
        HAS_TRANSCRIPT_AMPLIFICATION_CONSEQUENCE,
        HAS_INFRAME_INSERTION_CONSEQUENCE,
        HAS_INFRAME_DELETION_CONSEQUENCE,
        HAS_MISSENSE_VARIANT_CONSEQUENCE,
        HAS_PROTEIN_ALTERING_VARIANT_CONSEQUENCE,
        HAS_SPLICE_REGION_VARIANT_CONSEQUENCE,
        HAS_INCOMPLETE_TERMINAL_CODON_VARIANT_CONSEQUENCE,
        HAS_STOP_RETAINED_VARIANT_CONSEQUENCE,
        HAS_SYNONYMOUS_VARIANT_CONSEQUENCE,
        HAS_CODING_SEQUENCE_VARIANT_CONSEQUENCE,
        HAS_MATURE_MIRNA_VARIANT_CONSEQUENCE,
        HAS_5_PRIME_UTR_VARIANT_CONSEQUENCE,
        HAS_3_PRIME_UTR_VARIANT_CONSEQUENCE,
        HAS_NON_CODING_TRANSCRIPT_EXON_VARIANT_CONSEQUENCE,
        HAS_INTRON_VARIANT_CONSEQUENCE,
        HAS_NMD_TRANSCRIPT_VARIANT_CONSEQUENCE,
        HAS_NON_CODING_TRANSCRIPT_VARIANT_CONSEQUENCE,
        HAS_UPSTREAM_GENE_VARIANT_CONSEQUENCE,
        HAS_DOWNSTREAM_GENE_VARIANT_CONSEQUENCE,
        HAS_TFBS_ABLATION_CONSEQUENCE,
        HAS_TFBS_AMPLIFICATION_CONSEQUENCE,
        HAS_TF_BINDING_SITE_VARIANT_CONSEQUENCE,
        HAS_REGULATORY_REGION_ABLATION_CONSEQUENCE,
        HAS_REGULATORY_REGION_AMPLIFICATION_CONSEQUENCE,
        HAS_FEATURE_ELONGATION_CONSEQUENCE,
        HAS_REGULATORY_REGION_VARIANT_CONSEQUENCE,
        HAS_FEATURE_TRUNCATION_CONSEQUENCE,
        HAS_INTERGENIC_VARIANT_CONSEQUENCE
    }
    public void startDatabase(){
        log.log(Level.INFO, "Starting database ...");

        //create DB
        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }
    public void createIndexes() {
        Neo4j.createConstraint(graphDb, patientLabel, "PatientId");
        Neo4j.createConstraint(graphDb, sampleLabel, "SampleId");
        Neo4j.createConstraint(graphDb, snpLabel, "VariantId");
        Neo4j.createConstraint(graphDb, indelLabel, "VariantId");
        Neo4j.createConstraint(graphDb, symbolLabel, "SymbolId");
        Neo4j.createConstraint(graphDb, featureLabel, "FeatureId");
        Neo4j.createConstraint(graphDb, canonicalLabel, "FeatureId");
    }

    public void loadVCFFile(){
        log.log(Level.INFO, "Loading VCF into memory ...");

        int n = 0;

        Iterator<VariantContext> variantContextIterator = variantVcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            if (!variant.isFiltered() && variant.isVariant() && n < 2000){
                variants.add(variant);
                n++;
            }

        }
    }

    public void addPatientNodes(){
        log.log(Level.INFO, "Adding patient nodes ...");
        //TODO
    }
    public void addSampleNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample nodes ...");

        for (String sampleId : variantVcfFileReader.getFileHeader().getSampleNamesInOrder()){
            sampleNodeIds.put(sampleId, Neo4j.matchOrCreateUniqueNode(graphDb, sampleLabel, "SampleId", sampleId));
        }

    }
    public void extractVepAnnotations() throws InvalidPropertiesFormatException{
        log.log(Level.INFO, "Extracting functional annotations ...");
        String variantLookup;

        for (VariantContext variant : variants){

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotation vepAnnotation = new VEPAnnotation((String) variant.getAttribute("CSQ"));
                vepAnnotation.parseAnnotation();

                variantLookup = variant.getContig() + ":" + variant.getStart() + vepAnnotation.getAllele();

                //transcript, ensg, hgnc
                if (vepAnnotation.getFeature() != null && vepAnnotation.getGene() != null && vepAnnotation.getSymbol() != null) {

                    if(!vepAnnotations.containsKey(variantLookup)){
                        vepAnnotations.put(variantLookup, new HashSet<VEPAnnotation>());
                    }

                    vepAnnotations.get(variantLookup).add(vepAnnotation);

                }

            } catch (ClassCastException e) {

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")) {

                    VEPAnnotation vepAnnotation = new VEPAnnotation(annotation);
                    vepAnnotation.parseAnnotation();

                    variantLookup = variant.getContig() + ":" + variant.getStart() + vepAnnotation.getAllele();

                    //transcript, ensg, hgnc
                    if (vepAnnotation.getFeature() != null && vepAnnotation.getGene() != null && vepAnnotation.getSymbol() != null) {

                        if(!vepAnnotations.containsKey(variantLookup)){
                            vepAnnotations.put(variantLookup, new HashSet<VEPAnnotation>());
                        }

                        vepAnnotations.get(variantLookup).add(vepAnnotation);

                    }

                }

            }
        }
    }
    public void addVariantNodesAndGenotypeRelationships() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding variants and genotypes ...");

        GenomeVariant genomeVariant;

        //read VCF records
        for (VariantContext variant : variants) {
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                if (genotype.isNoCall() || genotype.isHomRef()) continue;
                if (genotype.getPloidy() != 2 || genotype.getAlleles().size() != 2) throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is not diploid");
                if (genotype.isMixed()) continue;

                if (genotype.isHomVar()) {

                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genotype.getSampleName(), genomeVariant, genotype.getGQ(), relTypes.HAS_HOM_VARIANT);

                } else if (genotype.isHetNonRef()) {

                    //het variant1
                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genotype.getSampleName(), genomeVariant, genotype.getGQ(), relTypes.HAS_HET_VARIANT);

                    //het variant2
                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genotype.getSampleName(), genomeVariant, genotype.getGQ(), relTypes.HAS_HET_VARIANT);

                } else if (genotype.isHet()) {

                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genotype.getSampleName(), genomeVariant, genotype.getGQ(), relTypes.HAS_HET_VARIANT);

                } else {
                    throw new InvalidPropertiesFormatException("Zygosity unknown: " + variant.toString());
                }

            }

        }

    }
    public void addPopulationFrequencies(){
        log.log(Level.INFO, "Adding population frequencies ...");
        //todo performance

        boolean foundPopulationFreq;
        GenomeVariant tempVariant;
        HashMap<String, Object> populationFrequencies = new HashMap<>();
        ArrayList<String> infoAtrributes;

        //loop over added variants
        for (Map.Entry<GenomeVariant, Long> genomeVariant : variantNodeIds.entrySet()){
            Iterator<VariantContext> variantContextIterator = oneKgP3VcfFileReader.query(genomeVariant.getKey().getContig(), genomeVariant.getKey().getPos() - 100, genomeVariant.getKey().getPos() + 100);

            foundPopulationFreq = false;

            //get 1kg variants
            loops:
            while (variantContextIterator.hasNext()) {
                VariantContext vcfVariant = variantContextIterator.next();

                for (int n = 0; n < vcfVariant.getAlternateAlleles().size(); ++n){

                    tempVariant = new GenomeVariant(vcfVariant.getContig(), vcfVariant.getStart(), vcfVariant.getReference().getBaseString(), vcfVariant.getAlternateAllele(n).getBaseString());
                    tempVariant.convertToMinimalRepresentation();

                    if (tempVariant.equals(genomeVariant.getKey())){

                        try{

                            //found variant
                            populationFrequencies.put("EAS_AF", Double.parseDouble((String) vcfVariant.getAttribute("EAS_AF")));
                            populationFrequencies.put("EUR_AF", Double.parseDouble((String) vcfVariant.getAttribute("EUR_AF")));
                            populationFrequencies.put("AFR_AF", Double.parseDouble((String) vcfVariant.getAttribute("AFR_AF")));
                            populationFrequencies.put("AMR_AF", Double.parseDouble((String) vcfVariant.getAttribute("AMR_AF")));
                            populationFrequencies.put("SAS_AF", Double.parseDouble((String) vcfVariant.getAttribute("SAS_AF")));

                        } catch (ClassCastException e) {

                            //found variant
                            infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("EAS_AF");
                            populationFrequencies.put("EAS_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                            infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("EUR_AF");
                            populationFrequencies.put("EUR_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                            infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("AFR_AF");
                            populationFrequencies.put("AFR_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                            infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("AMR_AF");
                            populationFrequencies.put("AMR_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                            infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("SAS_AF");
                            populationFrequencies.put("SAS_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                        }

                        Neo4j.addNodeProperty(graphDb, genomeVariant.getValue(), populationFrequencies);
                        populationFrequencies.clear();

                        foundPopulationFreq = true;
                        break loops;
                    }

                }

            }

            if(!foundPopulationFreq){
                log.log(Level.WARNING, genomeVariant.getKey().getConcatenatedVariant() + " not found in poly VCF");
            }

        }

    }
    public void addSymbolNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding symbols ...");

        for (Map.Entry<String, HashSet<VEPAnnotation>> variant : vepAnnotations.entrySet()){
            for (VEPAnnotation annotation : variant.getValue()){
                symbolNodeIds.put(annotation.getSymbol(), Neo4j.matchOrCreateUniqueNode(graphDb, symbolLabel, "SymbolId", annotation.getSymbol()));
            }
        }

    }
    public void addFeatureNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding features ...");
        ArrayList<Node> nodes;

        for (Map.Entry<String, HashSet<VEPAnnotation>> variant : vepAnnotations.entrySet()){
            for (VEPAnnotation annotation : variant.getValue()){

                nodes = Neo4j.getNodes(graphDb, featureLabel, "FeatureId", annotation.getFeature());
                if (nodes.size() == 0) {

                    try (Transaction tx = graphDb.beginTx()) {

                        Node node = graphDb.createNode();
                        node.addLabel(featureLabel);

                        node.setProperty("FeatureId", annotation.getFeature());
                        if (annotation.getFeatureType() != null) node.setProperty("FeatureType", annotation.getFeatureType());
                        if (annotation.getBiotype() != null) node.setProperty("Biotype", annotation.getBiotype());
                        if (annotation.getCanonical() != null && annotation.getCanonical().equals("YES")) node.addLabel(canonicalLabel);

                        featureNodeIds.put(annotation.getFeature(), node.getId());

                        tx.success();
                    }

                } else {
                    featureNodeIds.put(annotation.getFeature(), nodes.get(0).getId());
                }
            }
        }
    }
    public void addFunctionalAnnotationNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding functional annotations ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Long> variant : variantNodeIds.entrySet()){

            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())){
                if (getAnnotationsForVariant(variant.getKey(), graphDb).size() == 0){

                    if (!annotationNodeIds.containsKey(variant.getKey())) annotationNodeIds.put(variant.getKey(), new HashMap<String, Long>());

                    //loop over functional annotations for this variant
                    for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())) {

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(annotationLabel);

                            if(annotation.getExon() != null) {
                                String[] fields = annotation.getExon().split("/");
                                node.setProperty("Exon", fields[0]);
                                node.setProperty("TotalExons", fields[1]);
                            }
                            if(annotation.getHgvsCoding() != null) node.setProperty("HGVSc", annotation.getHgvsCoding());
                            if(annotation.getHgvsProtein() != null) node.setProperty("HGVSp", annotation.getHgvsProtein());
                            if(annotation.getIntron() != null) {
                                String[] fields = annotation.getIntron().split("/");
                                node.setProperty("Intron", fields[0]);
                                node.setProperty("TotalIntrons", fields[1]);
                            }
                            if(annotation.getStrand() != null) node.setProperty("Strand", annotation.getStrand());
                            if(annotation.getPolyphen() != null) node.setProperty("PolyPhen", annotation.getPolyphen());
                            if(annotation.getSift() != null) node.setProperty("SIFT", annotation.getSift());

                            annotationNodeIds.get(variant.getKey()).put(annotation.getFeature(), node.getId());

                            tx.success();
                        }

                    }
                }

            } else {
                log.log(Level.WARNING, "No annotation for: " + variant.getKey().getConcatenatedVariant());
            }
        }

    }
    public void addConsequenceRelationships() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Linking variants To functional annotations ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Long> variant : variantNodeIds.entrySet()){

            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())){

                HashMap<String, Object> properties = new HashMap<>();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())){

                    if (annotation.getConsequences().contains("transcript_ablation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TRANSCRIPT_ABLATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("splice_acceptor_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SPLICE_ACCEPTOR_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("splice_donor_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SPLICE_DONOR_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("stop_gained")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_STOP_GAINED_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("frameshift_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_FRAMESHIFT_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("stop_lost")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_STOP_LOST_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("start_lost")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_START_LOST_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("transcript_amplification")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TRANSCRIPT_AMPLIFICATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("inframe_insertion")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INFRAME_INSERTION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("inframe_deletion")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INFRAME_DELETION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("missense_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_MISSENSE_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("protein_altering_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_PROTEIN_ALTERING_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("splice_region_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SPLICE_REGION_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("incomplete_terminal_codon_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INCOMPLETE_TERMINAL_CODON_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("stop_retained_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_STOP_RETAINED_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("synonymous_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SYNONYMOUS_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("coding_sequence_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_CODING_SEQUENCE_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("mature_miRNA_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_MATURE_MIRNA_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("5_prime_UTR_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_5_PRIME_UTR_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("3_prime_UTR_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_3_PRIME_UTR_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("non_coding_transcript_exon_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_NON_CODING_TRANSCRIPT_EXON_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("intron_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INTRON_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("NMD_transcript_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_NMD_TRANSCRIPT_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("non_coding_transcript_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_NON_CODING_TRANSCRIPT_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("upstream_gene_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_UPSTREAM_GENE_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("downstream_gene_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_DOWNSTREAM_GENE_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("TFBS_ablation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TFBS_ABLATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("TFBS_amplification")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TFBS_AMPLIFICATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("TF_binding_site_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TF_BINDING_SITE_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("regulatory_region_ablation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_ABLATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("regulatory_region_amplification")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_AMPLIFICATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("feature_elongation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_FEATURE_ELONGATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("regulatory_region_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_VARIANT_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("feature_truncation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_FEATURE_TRUNCATION_CONSEQUENCE, properties);
                    if (annotation.getConsequences().contains("intergenic_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variant.getKey()), annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INTERGENIC_VARIANT_CONSEQUENCE, properties);

                }
            }

        }
    }
    public void addInFeatureRelationships() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Linking annotations to features ...");

        ArrayList<Node> nodes;

        //loop over variants
        for (Map.Entry<GenomeVariant, Long> variant : variantNodeIds.entrySet()){

            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())){

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())) {

                    //add in feature relationship
                    nodes = getSymbolsForFeature(annotation.getFeature(), graphDb);
                    if (nodes.size() == 0){

                        //link annotation to feature
                        try (Transaction tx = graphDb.beginTx()) {

                            Node node1 = graphDb.getNodeById(annotationNodeIds.get(variant.getKey()).get(annotation.getFeature()));
                            Node node2 = graphDb.getNodeById(featureNodeIds.get(annotation.getFeature()));

                            Relationship relationship = node1.createRelationshipTo(node2, relTypes.IN_FEATURE);
                            tx.success();

                        }

                    }

                }

            } else {
                log.log(Level.WARNING, "No annotation for: " + variant.getKey().getConcatenatedVariant());
            }
        }

    }
    public void addInSymbolRelationships() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Linking features to symbol ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Long> variant : variantNodeIds.entrySet()) {

            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())) {

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())) {

                    ArrayList<Node> nodes = getSymbolsForFeature(annotation.getFeature(), graphDb);

                    //check feature is not already associated with a symbol
                    if (nodes.size() == 0){
                        Neo4j.createRelationship(graphDb, featureNodeIds.get(annotation.getFeature()), symbolNodeIds.get(annotation.getSymbol()), relTypes.IN_SYMBOL, new HashMap<String, Object>());
                    }

                }
            }

        }

    }

    private static ArrayList<Node> getAnnotationsForVariant(GenomeVariant genomeVariant, final GraphDatabaseService graphDb){
        //todo convert to embedded
        ArrayList<Node> nodes = new ArrayList<>();

        for (Map<String, Object> result : Neo4j.runCypherQuery(graphDb,
                "MATCH (v:Variant {Variant:\"" + genomeVariant.getConcatenatedVariant() + "\"})-[]-(a:Annotation) return distinct a;"
        )){
            for (Map.Entry<String, Object> iter : result.entrySet()){
                nodes.add((Node) iter.getValue());
            }
        }

        return nodes;
    }
    private static ArrayList<Node> getSymbolsForFeature(String feature, final GraphDatabaseService graphDb){
        //todo convert to embedded

        ArrayList<Node> nodes = new ArrayList<>();

        for (Map<String, Object> result : Neo4j.runCypherQuery(graphDb,
                "MATCH (:Feature {FeatureId:\"" + feature + "\"})-[]-(s:Symbol) return distinct s;"
        )){
            for (Map.Entry<String, Object> iter : result.entrySet()){
                nodes.add((Node) iter.getValue());
            }
        }

        return nodes;
    }
    private static boolean doesGenotypeExist(String sampleId, String libraryId, String variantId, final GraphDatabaseService graphDb){

        //TODo convert to embedded function
        //todo allow snp/indel
        for (Map<String, Object> result : Neo4j.runCypherQuery(graphDb, "MATCH (:Sample {SampleId:\"" + sampleId + "\"})-[r]-(:Variant {VariantId:\"" + variantId + "\"}) where r.LibraryId = \"" + libraryId + "\" return count(r) as genotypes;")){
            if((Long) result.get("genotypes") > 0){
                return true;
            }
        }

        return false;
    }
    private void addVariantNodesAndGenotypeRelationshipsHelper(String sampleName, GenomeVariant genomeVariant, int genotypeQuality, RelationshipType relationshipType) throws InvalidPropertiesFormatException {

        if (!variantNodeIds.containsKey(genomeVariant)){

            if (genomeVariant.getRef().length() == 1 && genomeVariant.getAlt().length() == 1) {
                variantNodeIds.put(genomeVariant, Neo4j.matchOrCreateUniqueNode(graphDb, snpLabel, "VariantId", genomeVariant.getConcatenatedVariant()));
            } else if (genomeVariant.getRef().length() != 1 || genomeVariant.getAlt().length() != 1) {
                variantNodeIds.put(genomeVariant, Neo4j.matchOrCreateUniqueNode(graphDb, indelLabel, "VariantId", genomeVariant.getConcatenatedVariant()));
            } else {
                throw new InvalidPropertiesFormatException("Variant " + genomeVariant.getConcatenatedVariant() + " not SNP or Indel");
            }

        }

        if (!doesGenotypeExist(sampleName, libraryId, genomeVariant.getConcatenatedVariant(), graphDb)){
            HashMap<String, Object> properties = new HashMap<>();

            properties.put("GQ", genotypeQuality);
            properties.put("LibraryId", libraryId);

            Neo4j.createRelationship(graphDb, sampleNodeIds.get(sampleName), variantNodeIds.get(genomeVariant), relationshipType, properties);

        }

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}