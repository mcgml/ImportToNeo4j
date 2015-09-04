package nhs.genetics.cardiff;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 03/09/15.
 */
public class GeneMap2 {

    private static final Logger log = Logger.getLogger(GeneMap2.class.getName());
    private File morbidMapFilePath;
    private HashMap<String, HashSet<String>> geneMap2 = new HashMap<>();

    public GeneMap2(File morbidMapFilePath){
        this.morbidMapFilePath = morbidMapFilePath;
    }

    public void parseMorbidMap(){

        String line;

        try (BufferedReader reader = new BufferedReader(new FileReader(morbidMapFilePath))){

            //read file
            while ((line = reader.readLine()) != null) {

                String[] fields = line.split("\\|");

                if (fields.length < 12){
                    continue;
                }

                String[] geneFields = fields[5].split(",");
                String[] disorderFields = fields[11].split(";");

                for (String geneField : geneFields){
                    String geneFieldTrimmed = geneField.trim();

                    if (geneFieldTrimmed.equals("")) continue;

                    if (!geneMap2.containsKey(geneFieldTrimmed)) {
                        geneMap2.put(geneFieldTrimmed, new HashSet<String>());
                    }

                    for (String mim : disorderFields) {
                        String mimTrimmed = mim.trim();

                        if (mimTrimmed.equals("")) continue;

                        geneMap2.get(geneFieldTrimmed).add(mimTrimmed);

                    }

                }

            }

        } catch (IOException e){
            log.log(Level.SEVERE, "Could not read morbid map file: " + e.getMessage());
        }
    }

    public HashMap<String, HashSet<String>> getGeneMap2() {
        return geneMap2;
    }
}
