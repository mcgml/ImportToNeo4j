package nhs.genetics.cardiff;

import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.graphdb.ConstraintViolationException;
import org.neo4j.io.fs.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.InvalidPropertiesFormatException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    //todo update dbSNP

    private static final Logger log = Logger.getLogger(Main.class.getName());

    private static final String version = "1.0.1";
    private static boolean newDatabase = false, addAnnotations = false;

    public static void main(String[] args) throws InvalidPropertiesFormatException {

        if (args.length < 2 || args.length > 3) {
            System.err.println("ImportToNeo4j v" + version);
            System.err.println("Usage: <VCF> <db>");
            System.err.println("Options: -n New database, -a Annotated VCF");
            System.exit(1);
        }

        log.log(Level.INFO, "ImportToNeo4j v" + version);

        //update or overwrite, genotype or annotations?
        for (String arg : args){
            if (arg.equals("-n")){
                newDatabase = true;
            } else if (arg.equals("-a")){
                addAnnotations = true;
            }
        }

        if (newDatabase && addAnnotations){
            log.log(Level.SEVERE, "Cannot create new database and add annotations simultaneously. Check arguments.");
            System.exit(1);
        }

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
        VCFFileReader vcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));

        //create database object
        VariantDatabase variantDatabase = new VariantDatabase(vcfFileReader, new File(args[1]));
        variantDatabase.startDatabase();

        //add genotypes
        if (!addAnnotations){

            if (newDatabase) variantDatabase.createIndexes();

            try {
                variantDatabase.addSampleAndRunInfoNodes();
            } catch (ConstraintViolationException e){
                log.log(Level.SEVERE, "One or more analyses already exist in the database, check input.");
                System.exit(1);
            }

            variantDatabase.importVariants();
            variantDatabase.writeNewVariantsToVCF();

        } else {
            variantDatabase.importAnnotations();
        }

        variantDatabase.shutdownDatabase();
        vcfFileReader.close();

    }

}
