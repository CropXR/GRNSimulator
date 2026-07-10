
package ch.epfl.lis.gnw;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.apache.commons.logging.Log;

import ch.epfl.lis.imod.ImodNetwork;
import ch.epfl.lis.networks.Edge;
import ch.epfl.lis.networks.Node;
import ch.epfl.lis.networks.ios.ParseException;


/** Implements a Network creator that can create different types of network (for now random or modular and scale-free)
 */

public class NetworkCreator {
	
	
    private Logger log_ = Logger.getLogger(SubnetExtractor.class.getName());
    private NetworkType networkType_;
    private Integer numGenes_;
    private Integer numModules_;
    private Double sparsity_;
    private Double delta_in_;
    private Double delta_out_;
    private Double grouping_bias_;
    private String output_path_;
    
    private class MyEdge
    {
    	public MyEdge(Integer in, Integer out) 
    	{
    		Source = in;
    		Target = out;
    	}

		@Override
		public int hashCode() 
		{
			return Objects.hash(Source, Target);
		}
		
		@Override
		public boolean equals(Object obj)
		{
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			
			MyEdge other = (MyEdge) obj;
			return Objects.equals(Source, other.Source) && Objects.equals(Target, other.Target);
		}

		public Integer Source;
    	public Integer Target;
    };
    
    
    private ImodNetwork createRandomNetwork() throws ParseException
    {
		ImodNetwork network = new ImodNetwork();
		Random seed = new Random();
		HashSet<MyEdge> edges = new HashSet<MyEdge>();
        HashSet<Integer> vertices = new HashSet<>();

        while (vertices.size() < numGenes_) 
  		{
              double r = seed.nextDouble();
              int v, w;

              if( r < sparsity_)
              {
            	  v = vertices.size();
            	  vertices.add(v);
            	  // pick w randomly
            	  Double[] p = new Double[v + 1];
            	  Arrays.fill(p, 1.0);
               	  w = weightedRandomChoice(Arrays.asList(p), seed);
              }
              else
              {
            	  // pick v and w randomly
            	  v = vertices.size();
            	  Double[] p = new Double[v];
            	  Arrays.fill(p, 1.0);
               	  w = weightedRandomChoice(Arrays.asList(p), seed);
               	  v = weightedRandomChoice(Arrays.asList(p), seed);              
              }
              
              // add edge (w, v)
              MyEdge newEdge = new MyEdge(w, v);
              if ( !edges.contains(newEdge))
              {
            	  edges.add(newEdge);
              }
  		}
        
        GeneNetwork geneNetwork = createGeneNetwork(network, edges, vertices, seed);

        return geneNetwork; 
    }
    
    
    private ImodNetwork createModularScaleFreeNetwork() throws ParseException, IOException, URISyntaxException
    {
		ImodNetwork network = new ImodNetwork();
		network.setSigned(true);
		HashSet<MyEdge> edges = new HashSet<MyEdge>();
        HashSet<Integer> vertices = new HashSet<>(Arrays.asList(0, 1, 2)); 
		HashMap<Integer, HashSet<Integer>> modules = new HashMap<>();
		Random seed;
		URL url = GnwSettings.getInstance().getURL(output_path_ + "/" + network.getId() + "gene_modules.tsv");
		FileWriter fw = new FileWriter(new File(url.toURI()), false);
    	fw.write("Gene\tModule\n");

		if(GnwSettings.getInstance().getRandomSeed() == -1)
		{
			seed = new Random();
		}
		else
		{
			 seed = new Random(GnwSettings.getInstance().getRandomSeed()); 
		}
		
		if (numModules_ == 1)
		{
			modules.put(0, new HashSet<>(Arrays.asList(0, 1, 2)));
        	fw.write("0\t0\n");
        	fw.write("1\t0\n");
        	fw.write("2\t0\n");
		}
		else if (numModules_ == 2)
		{
			modules.put(0, new HashSet<>(Arrays.asList(0, 2)));
			modules.put(1, new HashSet<>(Arrays.asList(1)));
        	fw.write("0\t0\n");
        	fw.write("1\t1\n");
        	fw.write("2\t0\n");
		}
		else
		{
			modules.put(0, new HashSet<>(Arrays.asList(0)));
			modules.put(1, new HashSet<>(Arrays.asList(1)));	
			modules.put(2, new HashSet<>(Arrays.asList(2)));	
        	fw.write("0\t0\n");
        	fw.write("1\t1\n");
        	fw.write("2\t2\n");
		}
				
		//Start with small connected component
        List<Integer> edge_sources = new ArrayList<>(vertices);
        List<Integer> edge_targets = vertices.stream().map(i -> (i + 1) % vertices.size()).collect(Collectors.toList()); // edge targets
        
    	edges.add(new MyEdge(0, 1));
    	edges.add(new MyEdge(1, 2));
    	edges.add(new MyEdge(2, 0));

        HashMap<String, HashMap<Integer, Integer>> degrees = new HashMap<>();
        degrees.put("in", new HashMap<>());
        degrees.put("out", new HashMap<>());

        for (int i : vertices) 
        {
        	degrees.get("in").put(i, Collections.frequency(edge_targets, i)); // in-degree
        	degrees.get("out").put(i, Collections.frequency(edge_sources, i)); // out-degree
        }
        
        while (vertices.size() < numGenes_) 
		{
            double r = seed.nextDouble();
            Integer v, w;

            if( r < sparsity_)
            {
            	// gamma: add new node v with random class k_v
            	v = vertices.size();
            	vertices.add(v);
            	degrees.get("out").put(v, 0);
            	degrees.get("in").put(v, 0);
            	Integer k_v = seed.nextInt(numModules_);
            	modules.computeIfAbsent(k_v, k -> new HashSet<>()).add(v);

            	// pick w by out degree biased towards type k_v
            	ArrayList<Double> p = new ArrayList<>();
            	for (int n : vertices) 
            	{
            		p.add((modules.get(k_v).contains(n) ? grouping_bias_ : 1) * (degrees.get("out").get(n) + delta_out_));
            	}
            	w = weightedRandomChoice(p, seed);
            	fw.write(v.toString() + "\t" + k_v.toString() + "\n");
            }
            else 
            {
                // beta: pick w by out degree with bias delta_out and v as in alpha regime
            	ArrayList<Double> p = new ArrayList<>();
                for (int n : vertices) 
                {
                    p.add(degrees.get("out").get(n) + delta_out_);
                }
                w = weightedRandomChoice(p, seed);
                int k_w = modules.entrySet().stream().filter(entry -> entry.getValue().contains(w)).findFirst().get().getKey();
                ArrayList<Double> q = new ArrayList<>();
                for (int n : vertices) 
                {
                    q.add((modules.get(k_w).contains(n) ? grouping_bias_ : 1) * (degrees.get("in").get(n) + delta_in_));
                }
                v = weightedRandomChoice(q, seed);
            }
            
            // add edge (w, v)
            MyEdge newEdge = new MyEdge(w, v);
            if ( !edges.contains(newEdge))
            {
            	edge_sources.add(w);
            	edge_targets.add(v);
            	degrees.get("out").put(w, degrees.get("out").get(w) + 1);
            	degrees.get("in").put(v, degrees.get("in").get(v) + 1);
            	edges.add(newEdge);
            }
 		}
             
        GeneNetwork geneNetwork = createGeneNetwork(network, edges, vertices, seed);
        fw.close();
		return geneNetwork;
    }


	private GeneNetwork createGeneNetwork(ImodNetwork network, HashSet<MyEdge> edges, HashSet<Integer> vertices,
			Random seed) throws ParseException
	{
		for ( int v : vertices)
        {
        	network.addNode(new Node("Node_" + v));
        }
        for ( MyEdge edge : edges)
        {
        	int source = edge.Source;
        	int target = edge.Target;
            double r = seed.nextDouble();
            if( r < 0.5 )
            {
            	network.addEdge(new Edge(network.getNode(source), network.getNode(target), "+" ));
            }
            else
            {
            	network.addEdge(new Edge(network.getNode(source), network.getNode(target), "-" ));
            }
        }
        
       
        network.setId("");
        GeneNetwork geneNetwork = new GeneNetwork(network);
		return geneNetwork;
	}
    
    private int weightedRandomChoice(List<Double> p, Random seed)
    {
    	Double total = p.stream().mapToDouble(Double::doubleValue).sum();
    	Double value = total * seed.nextDouble();
    	Double curr = 0.;
    	for ( int indx = 0; indx < p.size(); ++indx)
    	{
    		curr += p.get(indx);
    		if( curr > value )
    		{
    			return indx;
    		}	
    	}
    	
		return  p.size() - 1;
	}

	public enum NetworkType 
    {
  	  RANDOM,
      MODULAR_SCALEFREE,	
  	  SMALLWORLD,
  	  SCALEFREE,
  	  SMALLWORLD_SCALEFREE
  	}
    
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public NetworkCreator(NetworkType type, Integer numGenes, Integer numModules, Double sparsity, Double delta_in, Double delta_out, Double grouping_bias, String outputPath ) 
	{
		networkType_ = type;
		numGenes_ = numGenes;
		numModules_ = numModules;
		sparsity_ = sparsity;
		delta_in_ = delta_in;
		delta_out_ = delta_out;
		grouping_bias_ = grouping_bias;
		output_path_ = outputPath;
	}
	
	public ImodNetwork createNetwork() throws ParseException, IllegalArgumentException, IOException, URISyntaxException
	{
		ImodNetwork network = new ImodNetwork();
		switch (networkType_)
		{
		
		case RANDOM: 
			//create a random network 
			network = createRandomNetwork();
			break;
		case MODULAR_SCALEFREE: 
			//Create modular scale-free network following https://www.biorxiv.org/content/10.1101/2024.07.04.602130v2		
			network = createModularScaleFreeNetwork();	
			break;
			
		case SMALLWORLD:
			
		case SCALEFREE:
			
		case SMALLWORLD_SCALEFREE:
			
		 default:
			 log_.log(Level.SEVERE, "NetworkType doesn't exist");
			 throw new IllegalArgumentException("Unsupported network type");
		}
		
		return network;
	}
}
	