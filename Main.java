package nhs.genetics.cardiff;

import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.io.fs.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.InvalidPropertiesFormatException;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final double version = 0.1;

    public static void main(String[] args) throws InvalidPropertiesFormatException {

        if (args.length != 5) {
            System.err.println("Usage: <GenotypeVCFFile> <AnnotatedVCFFile> <GeneMap2File> <Neo4jDBPath> <Overwrite>");
            System.exit(1);
        }

        log.log(Level.INFO, "ImportToNeo4j v" + version);
        log.log(Level.INFO, "Importing " + args[0] + " to " + args[3]);

        //create VCF file parser
        VCFFileReader variantVcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));
        VCFFileReader annotationVcfFileReader = new VCFFileReader(new File(args[1]), new File(args[1] + ".idx"));

        //create OMIM file parser
        GeneMap2 geneMap2 = new GeneMap2(new File(args[2]));
        geneMap2.parseMorbidMap();

        //create database object
        VariantDatabase variantDatabase = new VariantDatabase(variantVcfFileReader, annotationVcfFileReader, new File(args[3]), geneMap2.getGeneMap2());

        //update or overwrite
        boolean overwriteDB = Boolean.parseBoolean(args[4]);

        //delete existing DB
        if (overwriteDB) {
            try{
                FileUtils.deleteRecursively(new File(args[3]));
            } catch (IOException e){
                log.log(Level.SEVERE, "Could not delete database: " + e.getMessage());
                System.exit(1);
            }
        }

        //create new DB
        variantDatabase.startDatabase();
        if (overwriteDB) variantDatabase.createIndexes();

        variantDatabase.loadVCFFiles();

        variantDatabase.addSampleAndRunInfoNodes();
        variantDatabase.addVariantNodesAndGenotypeRelationships();
        variantDatabase.addAnnotations();

        variantDatabase.shutdownDatabase();

        variantVcfFileReader.close();
        annotationVcfFileReader.close();

    }

}
