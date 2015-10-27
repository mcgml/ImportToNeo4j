package nhs.genetics.cardiff;

import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.io.fs.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.InvalidPropertiesFormatException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final double version = 0.3;
    private static boolean newDatabase = false;

    public static void main(String[] args) throws InvalidPropertiesFormatException {

        if (args.length < 2 || args.length > 3) {
            System.err.println("Usage: <VCF> <db>");
            System.err.println("Options: -n new database");
            System.exit(1);
        }

        log.log(Level.INFO, "ImportToNeo4j v" + version);

        //update or overwrite
        for (String arg : args){
            if (args.equals("-n")){
                newDatabase = true;
            }
        }

        log.log(Level.INFO, "New database: " + newDatabase);

        if (newDatabase) {
            log.log(Level.INFO, "Deleting existing database");
            try{
                FileUtils.deleteRecursively(new File(args[1]));
            } catch (IOException e){
                log.log(Level.SEVERE, "Could not delete database: " + e.getMessage());
                System.exit(1);
            }
        }

        log.log(Level.INFO, "Importing " + args[0] + " to " + args[1]);

        //create VCF file parser
        VCFFileReader variantVcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));

        //create database object
        VariantDatabase variantDatabase = new VariantDatabase(variantVcfFileReader, new File(args[1]));

        //create new DB
        variantDatabase.startDatabase();
        if (newDatabase) variantDatabase.createIndexes();

        variantDatabase.loadVCF();
        if (newDatabase) variantDatabase.addUsers();

        variantDatabase.addSampleAndRunInfoNodes();
        variantDatabase.addVariantNodesAndGenotypeRelationships();
        variantDatabase.addAnnotations();
        if (newDatabase) variantDatabase.addGenePanels();

        variantDatabase.shutdownDatabase();

        variantVcfFileReader.close();

    }

}
