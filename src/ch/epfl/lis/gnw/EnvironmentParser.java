package ch.epfl.lis.gnw;

import javax.xml.parsers.DocumentBuilderFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import ch.epfl.lis.imod.ImodNetwork;
import ch.epfl.lis.networks.Edge;

public class EnvironmentParser 
{
    private static Logger log_ = Logger.getLogger(GeneNetwork.class.getName());

	 public static void parse(String filename, ImodNetwork grn) 
	 {
		 try 
		 {
			 DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
	         DocumentBuilder builder = factory.newDocumentBuilder();
	         Document document = builder.parse(filename);
		     
		     NodeList envNodes = document.getElementsByTagName("Environment");

		     for (int i = 0; i < envNodes.getLength(); i++) 
		     {
		            Element envElem = (Element) envNodes.item(i);

		            // Extract attributes
		            String label = envElem.getAttribute("label");
		            String type = envElem.getAttribute("type");
		                       
		            //create the environment node
		            GrnNode node = null;
		            if( grn instanceof GeneNetwork)
		            {
		            	node = GrnNodeFactory.createNode(type, (GeneNetwork)grn);
			            parseParameters(envElem, node);	
		            }
		            else
		            {
		            	node = new HillGene((GeneNetwork)grn);
		            }
		            
		            node.setLabel(label);	            
	            	parseTargets(grn, envElem, node);		            
		            grn.addNode(node);
		     }
	     } 
		 catch (Exception e) {
	            e.printStackTrace();
	     }
	}

	 private static void parseTargets(ImodNetwork grn, Element envElem, GrnNode node) 
	 {
		ArrayList<GrnNode> targets = new ArrayList<GrnNode>();
		NodeList targetsNodes = envElem.getElementsByTagName("targets");
		if (targetsNodes.getLength() > 0) 
		{
			 String targetsText = targetsNodes.item(0).getTextContent().trim();
		     String[] targetsStr = targetsText.split("\\s+");
		     for (String t : targetsStr) 
		     {
		         if (!t.isEmpty())
		         {
		        	 targets.add((GrnNode)grn.getNode(Integer.parseInt(t)));
		         }
		     }
		}
		else if (envElem.getAttribute("numTargets") != null)
		{
			//choose numTargets random target
			int numTargets = Integer.parseInt(envElem.getAttribute("numTargets"));
			 
			if (numTargets > grn.getSize())
			{
				log_.log(Level.WARNING, "numTargets is larger than the total number of genes for Environment "  
						+ node.getLabel());
				numTargets = grn.getSize();
			}
			
			List<Integer> numbers = new ArrayList<>();
			for (int indx = 0; indx < grn.getSize(); indx++) {
				numbers.add(indx);
			}
		
			Collections.shuffle(numbers);
			List<Integer> targetIndices = numbers.subList(0, numTargets);
		
			//set these as target genes
			for(Integer targetIndx : targetIndices)
			{
				targets.add( (GrnNode)grn.getNode(targetIndx));
			}
		}
		else
		{				
			log_.log(Level.SEVERE, "Either specify numTargets or a target list for Environment "  
				+ node.getLabel());
			
		}
		
		for ( GrnNode target : targets)
		{
			grn.addEdge(new Edge(node, target, Edge.ENHANCER));
		}
	 }

	 private static void parseParameters(Element envElem, GrnNode node) 
	 {
		//Parse listOfParameters
		ArrayList<String> parameterNames = new ArrayList<String>(); // the names of the parameters
		ArrayList<Double> parameterValues = new ArrayList<Double>(); // the values of the parameters
		NodeList listOfParamsNodes = envElem.getElementsByTagName("listOfParameters");
		if (listOfParamsNodes.getLength() > 0)
		{
		    Element listOfParams = (Element) listOfParamsNodes.item(0);
		    NodeList paramNodes = listOfParams.getElementsByTagName("parameter");

		    for (int j = 0; j < paramNodes.getLength(); j++) 
		    {
		        Element param = (Element) paramNodes.item(j);

		        String id = param.getAttribute("id");
		        double value = Double.parseDouble(param.getAttribute("value"));
		        parameterNames.add(id);
		        parameterValues.add(value);
		    }              
		}
		
		node.initialization(parameterNames, parameterValues, new ArrayList<GrnNode>());
	 }
}
