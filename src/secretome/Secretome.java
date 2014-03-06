/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package secretome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import org.json.simple.parser.ContainerFactory;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

/**
 *
 * @author ethering
 */
public class Secretome
{

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, org.json.simple.parser.ParseException
    {
        File outfile = new File("/Users/ethering/projects/charala/go_found_not_found.txt");

        File secretome = new File("/Users/ethering/projects/charala/Charala_secretome_list.txt");
        String jsonFile = "/Users/ethering/git_repositories/h_pseu_analysis/go_analysis/terms_and_descriptions.json";
        Secretome jp = new Secretome();
        jp.createGoCounts(outfile, secretome, jsonFile);
    }

    public void createGoCounts(File outfile, File SecretomeList, String jsonTerms) throws IOException, ParseException
    {

        FileWriter fw = new FileWriter(outfile.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write("go_term\tdescription\tsecretome_genes\tnon_secretome_genes\tprob of combination\texpected\tpopulation" + System.getProperty("line.separator"));
        HashMap<String, ArrayList<String>> goMappings = new HashMap<>();
        //read in the protein ids for the secretome and add them to an array list
        Scanner s = new Scanner(SecretomeList);
        ArrayList<String> secretomeList = new ArrayList<>();
        while (s.hasNext())
        {
            secretomeList.add(s.next());

        }
        s.close();

        //open the json file and read into a string
        String jsonText = readFile(jsonTerms, Charset.defaultCharset());
        JSONParser parser = new JSONParser();
        ContainerFactory containerFactory = new ContainerFactory()
        {
            public List creatArrayContainer()
            {
                return new LinkedList();
            }

            public Map createObjectContainer()
            {
                return new LinkedHashMap();
            }
        };

        //create a hash of key (protein id),value (go terms) pairs from the json
        Map json = (Map) parser.parse(jsonText, containerFactory);
        Iterator iter = json.entrySet().iterator();
        //iterate over them
        while (iter.hasNext())
        {
            Map.Entry entry = (Map.Entry) iter.next();
            //get the gene id
            String gene = (String) entry.getKey();
            //get the go terms (GO:0000123, some_process)
            LinkedHashMap goTermsMap = (LinkedHashMap) entry.getValue();
            //System.out.println(goTermsMap);

            //iterate over the go ids only (ignore the process names)
            Set<String> goTerms = goTermsMap.keySet();

            for (String goTerm : goTerms)
            {
                String desc = (String) goTermsMap.get(goTerm);

                goTerm = goTerm.concat("\t" + desc);
                //System.out.println(goTerm);
                //if the go id is already in our hashmap, get the array list of gene ids and add the current gene id
                if (goMappings.containsKey(goTerm))
                {
                    ArrayList<String> al = goMappings.get(goTerm);
                    al.add(gene);
                    goMappings.put(goTerm, al);

                }
                // if it's not, create a new hashmap entry with the go id as the key and an arryaylist of gene ids as the value
                else
                {
                    ArrayList<String> al = new ArrayList<>();
                    al.add(gene);
                    goMappings.put(goTerm, al);
                }
            }
        }
        double noInSecretome = 0;
        double noNotInSecretome = 0;
        //iterate over the goMappings
        Iterator it = goMappings.entrySet().iterator();
        while (it.hasNext())
        {
            Map.Entry pairs = (Map.Entry) it.next();
            //get the list of protein ids..
            ArrayList<String> genesInGo = (ArrayList<String>) pairs.getValue();
            //System.out.println(goId);
            //..and see if they are in the secretome
            for (String geneInGo : genesInGo)
            {
                if (secretomeList.contains(geneInGo))
                {
                    noInSecretome++;
                }
                else
                {
                    noNotInSecretome++;
                }
            }
        }
        double probOfOneSecretomeGene = noInSecretome / (noInSecretome + noNotInSecretome);
        double probOfOneNonSecretomeGene = noNotInSecretome / (noInSecretome + noNotInSecretome);
        Iterator it2 = goMappings.entrySet().iterator();
        while (it2.hasNext())
        {
            double inSecretome = 0;
            double notInSecretome = 0;
            Map.Entry pairs = (Map.Entry) it2.next();
            String goId = (String) pairs.getKey();
            //get the list of protein ids..
            ArrayList<String> genesInGo = (ArrayList<String>) pairs.getValue();
            //System.out.println(goId);
            //..and see if they are in the secretome
            for (String geneInGo : genesInGo)
            {
                if (secretomeList.contains(geneInGo))
                {
                    inSecretome++;

                }
                else
                {
                    notInSecretome++;

                }
            }

            double totalPop = inSecretome + notInSecretome;

            double probOfCombination = inSecretome / totalPop;
            double expected = Math.pow(probOfOneSecretomeGene, inSecretome) * Math.pow(probOfOneNonSecretomeGene, notInSecretome);
            bw.write(goId + "\t" + inSecretome + "\t" + notInSecretome + "\t" + probOfCombination + "\t" + String.format(Locale.ENGLISH, "%.6f", expected) + "\t" + totalPop + System.getProperty("line.separator"));

        }
        bw.close();
        System.out.println("In secretome " + noInSecretome);
        System.out.println("Not in secretome " + noNotInSecretome);

    }

    static String readFile(String path, Charset encoding)
            throws IOException
    {
        byte[] encoded = Files.readAllBytes(Paths.get(path));
        return encoding.decode(ByteBuffer.wrap(encoded)).toString();
    }
}
