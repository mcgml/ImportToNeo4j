package nhs.genetics.cardiff;

import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 17/06/2015.
 */
public class VEPAnnotation {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    private String record, allele, gene, feature, featureType, cDNAPosition, cdsPosition, proteinPosition, aminoAcids, codons, existingVariation, aaMaf, eaMaf, exon, intron, distance, strand, symbol, sift, polyphen, gmaf, hgvsCoding, hgvsProtein, afrMaf, amrMaf, asnMaf, eurMaf;
    private HashSet<String> consequences = new HashSet<>();
    private HashSet<String> clinSigs = new HashSet<>();

    public VEPAnnotation(String record) {
        this.record = record;
    }

    public void parseAnnotation() {
        String[] fields = record.split("\\|");

        try {
            if (!fields[0].equals("")) {
                this.allele = fields[0];
            }

            if (!fields[1].equals("")) {
                this.gene = fields[1];
            }

            if (!fields[2].equals("")) {
                this.feature = fields[2];
            }

            if (!fields[3].equals("")) {
                this.featureType = fields[3];
            }

            if (!fields[4].equals("")) {
                for (String consequence : fields[4].split("&")){
                    this.consequences.add(consequence);
                }
            }

            if (!fields[5].equals("")) {
                this.cDNAPosition = fields[5];
            }

            if (!fields[6].equals("")) {
                this.cdsPosition = fields[6];
            }

            if (!fields[7].equals("")) {
                this.proteinPosition = fields[7];
            }

            if (!fields[8].equals("")) {
                this.aminoAcids = fields[8];
            }

            if (!fields[9].equals("")) {
                this.codons = fields[9];
            }

            if (!fields[10].equals("")) {
                this.existingVariation = fields[10];
            }

            if (!fields[11].equals("")) {
                this.aaMaf = fields[11];
            }

            if (!fields[12].equals("")) {
                this.eaMaf = fields[12];
            }

            if (!fields[13].equals("")) {
                this.exon = fields[13];
            }

            if (!fields[14].equals("")) {
                this.intron = fields[14];
            }

            if (!fields[15].equals("")) {
                this.distance = fields[15];
            }

            if (!fields[16].equals("")) {
                this.strand = fields[16];
            }

            if (!fields[17].equals("")) {
                for (String clinSig : fields[17].split("&")){
                    this.clinSigs.add(clinSig);
                }
            }

            if (fields[19].equals("HGNC") && !fields[18].equals("")) {
                this.symbol = fields[18];
            }

            if (!fields[20].equals("")) {
                String[] subFields = fields[20].split("\\(");
                this.sift = subFields[0];
            }

            if (!fields[21].equals("")) {
                String[] subFields = fields[21].split("\\(");
                this.polyphen = subFields[0];
            }

            if (!fields[22].equals("")) {
                this.gmaf = fields[22];
            }

            if (!fields[23].equals("")) {
                String[] subFields = fields[23].split(":");
                this.hgvsCoding = subFields[1];
            }

            if (!fields[24].equals("")) {
                String[] subFields = fields[24].split(":");

                if (subFields[1].contains("(") && subFields[1].contains(")")){
                    subFields = subFields[1].split("\\(");
                    subFields = subFields[1].split("\\)");
                    if (subFields[0].equals("p.%3D")) this.hgvsProtein = "p.="; else this.hgvsProtein = subFields[0];
                } else {
                    if (subFields[1].equals("p.%3D")) this.hgvsProtein = "p.="; else this.hgvsProtein = subFields[1];
                }
            }

            if (!fields[25].equals("")) {
                this.afrMaf = fields[25];
            }

            if (!fields[26].equals("")) {
                this.amrMaf = fields[26];
            }

            if (!fields[27].equals("")) {
                this.asnMaf = fields[27];
            }

            if (!fields[28].equals("")) {
                this.eurMaf = fields[28];
            }

        } catch (ArrayIndexOutOfBoundsException e) {
            log.log(Level.FINE, e.getMessage());
        }

    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VEPAnnotation that = (VEPAnnotation) o;

        if (allele != null ? !allele.equals(that.allele) : that.allele != null) return false;
        if (gene != null ? !gene.equals(that.gene) : that.gene != null) return false;
        if (feature != null ? !feature.equals(that.feature) : that.feature != null) return false;
        if (featureType != null ? !featureType.equals(that.featureType) : that.featureType != null) return false;
        if (cDNAPosition != null ? !cDNAPosition.equals(that.cDNAPosition) : that.cDNAPosition != null) return false;
        if (cdsPosition != null ? !cdsPosition.equals(that.cdsPosition) : that.cdsPosition != null) return false;
        if (proteinPosition != null ? !proteinPosition.equals(that.proteinPosition) : that.proteinPosition != null)
            return false;
        if (aminoAcids != null ? !aminoAcids.equals(that.aminoAcids) : that.aminoAcids != null) return false;
        if (codons != null ? !codons.equals(that.codons) : that.codons != null) return false;
        if (existingVariation != null ? !existingVariation.equals(that.existingVariation) : that.existingVariation != null)
            return false;
        if (aaMaf != null ? !aaMaf.equals(that.aaMaf) : that.aaMaf != null) return false;
        if (eaMaf != null ? !eaMaf.equals(that.eaMaf) : that.eaMaf != null) return false;
        if (exon != null ? !exon.equals(that.exon) : that.exon != null) return false;
        if (intron != null ? !intron.equals(that.intron) : that.intron != null) return false;
        if (distance != null ? !distance.equals(that.distance) : that.distance != null) return false;
        if (strand != null ? !strand.equals(that.strand) : that.strand != null) return false;
        if (symbol != null ? !symbol.equals(that.symbol) : that.symbol != null) return false;
        if (sift != null ? !sift.equals(that.sift) : that.sift != null) return false;
        if (polyphen != null ? !polyphen.equals(that.polyphen) : that.polyphen != null) return false;
        if (gmaf != null ? !gmaf.equals(that.gmaf) : that.gmaf != null) return false;
        if (hgvsCoding != null ? !hgvsCoding.equals(that.hgvsCoding) : that.hgvsCoding != null) return false;
        if (hgvsProtein != null ? !hgvsProtein.equals(that.hgvsProtein) : that.hgvsProtein != null) return false;
        if (afrMaf != null ? !afrMaf.equals(that.afrMaf) : that.afrMaf != null) return false;
        if (amrMaf != null ? !amrMaf.equals(that.amrMaf) : that.amrMaf != null) return false;
        if (asnMaf != null ? !asnMaf.equals(that.asnMaf) : that.asnMaf != null) return false;
        if (eurMaf != null ? !eurMaf.equals(that.eurMaf) : that.eurMaf != null) return false;
        if (consequences != null ? !consequences.equals(that.consequences) : that.consequences != null) return false;
        return !(clinSigs != null ? !clinSigs.equals(that.clinSigs) : that.clinSigs != null);

    }

    @Override
    public int hashCode() {
        int result = allele != null ? allele.hashCode() : 0;
        result = 31 * result + (gene != null ? gene.hashCode() : 0);
        result = 31 * result + (feature != null ? feature.hashCode() : 0);
        result = 31 * result + (featureType != null ? featureType.hashCode() : 0);
        result = 31 * result + (cDNAPosition != null ? cDNAPosition.hashCode() : 0);
        result = 31 * result + (cdsPosition != null ? cdsPosition.hashCode() : 0);
        result = 31 * result + (proteinPosition != null ? proteinPosition.hashCode() : 0);
        result = 31 * result + (aminoAcids != null ? aminoAcids.hashCode() : 0);
        result = 31 * result + (codons != null ? codons.hashCode() : 0);
        result = 31 * result + (existingVariation != null ? existingVariation.hashCode() : 0);
        result = 31 * result + (aaMaf != null ? aaMaf.hashCode() : 0);
        result = 31 * result + (eaMaf != null ? eaMaf.hashCode() : 0);
        result = 31 * result + (exon != null ? exon.hashCode() : 0);
        result = 31 * result + (intron != null ? intron.hashCode() : 0);
        result = 31 * result + (distance != null ? distance.hashCode() : 0);
        result = 31 * result + (strand != null ? strand.hashCode() : 0);
        result = 31 * result + (symbol != null ? symbol.hashCode() : 0);
        result = 31 * result + (sift != null ? sift.hashCode() : 0);
        result = 31 * result + (polyphen != null ? polyphen.hashCode() : 0);
        result = 31 * result + (gmaf != null ? gmaf.hashCode() : 0);
        result = 31 * result + (hgvsCoding != null ? hgvsCoding.hashCode() : 0);
        result = 31 * result + (hgvsProtein != null ? hgvsProtein.hashCode() : 0);
        result = 31 * result + (afrMaf != null ? afrMaf.hashCode() : 0);
        result = 31 * result + (amrMaf != null ? amrMaf.hashCode() : 0);
        result = 31 * result + (asnMaf != null ? asnMaf.hashCode() : 0);
        result = 31 * result + (eurMaf != null ? eurMaf.hashCode() : 0);
        result = 31 * result + (consequences != null ? consequences.hashCode() : 0);
        result = 31 * result + (clinSigs != null ? clinSigs.hashCode() : 0);
        return result;
    }

    public String getAllele() {
        return allele;
    }
    public String getGene() {
        return gene;
    }
    public String getFeature() {
        return feature;
    }
    public String getFeatureType() {
        return featureType;
    }
    public String getcDNAPosition() {
        return cDNAPosition;
    }
    public String getCdsPosition() {
        return cdsPosition;
    }
    public String getProteinPosition() {
        return proteinPosition;
    }
    public String getAminoAcids() {
        return aminoAcids;
    }
    public String getCodons() {
        return codons;
    }
    public String getExistingVariation() {
        return existingVariation;
    }
    public String getAaMaf() {
        return aaMaf;
    }
    public String getEaMaf() {
        return eaMaf;
    }
    public String getExon() {
        return exon;
    }
    public String getIntron() {
        return intron;
    }
    public String getDistance() {
        return distance;
    }
    public String getStrand() {
        return strand;
    }
    public String getSymbol() {
        return symbol;
    }
    public String getSift() {
        return sift;
    }
    public String getPolyphen() {
        return polyphen;
    }
    public String getGmaf() {
        return gmaf;
    }
    public String getHgvsCoding() {
        return hgvsCoding;
    }
    public String getHgvsProtein() {
        return hgvsProtein;
    }
    public String getAfrMaf() {
        return afrMaf;
    }
    public String getAmrMaf() {
        return amrMaf;
    }
    public String getAsnMaf() {
        return asnMaf;
    }
    public String getEurMaf() {
        return eurMaf;
    }
    public HashSet<String> getConsequences() {
        return consequences;
    }
    public HashSet<String> getClinSigs() {
        return clinSigs;
    }

}