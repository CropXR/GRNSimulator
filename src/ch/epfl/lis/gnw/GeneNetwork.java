/*
Copyright (c) 2008-2010 Daniel Marbach & Thomas Schaffter

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper(s) listed
on http://gnw.sourceforge.net.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

package ch.epfl.lis.gnw;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.jet.random.Uniform;
import ch.epfl.lis.networks.Edge;
import ch.epfl.lis.networks.Node;
import ch.epfl.lis.networks.Structure;
import ch.epfl.lis.imod.ImodNetwork;

import jigcell.sbml2.Compartment;
import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.SBMLLevel2Document;
import jigcell.sbml2.Species;
import jigcell.sbml2.ModifierSpeciesReference;


/** This class represents a gene network.
 * 
 * It extends ImodNetwork, thus it contains all the
 * genes and the edges. Furthermore, it implements the state (gene and protein expression
 * levels) and functions to compute the production rates and to load and save the network
 * in SBML format.
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 * @author Daniel Marbach (firstname.name@gmail.com)
 * 
 */ 
public class GeneNetwork extends ImodNetwork { 

	/** Current gene expression levels */
	private DoubleMatrix1D x_ = null;
	/** Current protein expression levels */
	private DoubleMatrix1D y_ = null;
	/** The number of regulators in the network (must be initialized using markTfs()) */
	private int numTfs_ = -1;
	
	private double time_ = 0;
    /** Logger for this class */
    private static Logger log_ = Logger.getLogger(GeneNetwork.class.getName());
	/** SBML file format (see load()) */
	public static final int SBML = 4;
		
	
	// ============================================================================
	// INITIALIZATION
	
	/**
	 * Default constructor
	 */
	public GeneNetwork() {
		super();
	}

	
	// ----------------------------------------------------------------------------

	/**
	 * Constructor from a network (copies the structure of the given network)
	 */
	public GeneNetwork(Structure network) 
	{

		// initialize with members of network, except for the nodes and the edges
		super(network);
		nodes_.clear();
		edges_.clear();
		
		for (int i=0; i < network.getSize(); i++) 
		{	
			HillGene g = new HillGene(this);
			g.setLabel(network.getNode(i).getLabel());
			nodes_.add(g);
		}
			
		Node source, target;
		Edge edge;
		for (int i=0; i < network.getNumEdges(); i++) {
			edge = network.getEdge(i);
			// Get the nodes in this.nodes_ and *not* in network.nodes_
			source = getNode(edge.getSource().getLabel());
			target = getNode(edge.getTarget().getLabel());
			edges_.add( new Edge(source, target, edge.getType()) ); //, edge.getLabel()) );
		}
		setSize(getSize());
	}
	

	// ============================================================================
	// PUBLIC METHODS
	
	/**
	 * Set network size, reserve space in x_ and y_ (does not affect members of the superclass)
	 * @param numGenes Size of the network
	 */
	public void setSize(int numGenes) {
		
		x_ = new DenseDoubleMatrix1D(numGenes);
		x_.assign(-1);
		if (GnwSettings.getInstance().getModelTranslation()) {
			y_ = new DenseDoubleMatrix1D(numGenes);
			y_.assign(-1);
		} else
			y_ = x_;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Based on W_, determine for each gene which are its inputs
	 */
	public void initializeInputWiring() {
		
		int size = getSize();
		
		for (int i=0; i < size; i++) {
			ArrayList<GrnNode> inputs = getInputs(nodes_.get(i));
			if( nodes_.get(i) instanceof Gene)
			{
				((Gene)nodes_.get(i)).setInputGenes(inputs);
			}
		}
	}

		
	// ----------------------------------------------------------------------------
	
	/**
	 * Random initialization of the dynamical model (based on the fixed topology defined by W_)
	 */
	public void randomInitialization() {
		
		log_.log(Level.INFO, "\nRandom initialization of " + id_);

		// set the inputs for the genes based on W_
		initializeInputWiring();
		int size = nodes_.size();

		Gene gene = null;
		for (int i=0; i < size; i++) 
		{
			if( nodes_.get(i) instanceof Gene)
			{
				gene = (Gene)nodes_.get(i);
				gene.randomInitialization();
				log_.log(Level.INFO, id_ + ": " + gene.getLabel() + ": " + gene.toString());
			}
		}
		// If the network was unsigned before the initialization, we have to set the types
		// of the edges (enhancing / inhibiting) accordingly. If the network was signed,
		// this just sets the same values that are already there.
		setEdgeTypesAccordingToDynamicalModel();
		signed_ = true;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Compute the rate of change of all state variables x (and y if translation is modelled).
	 * The result is returned concatenated into one array [dxdt dydt].
	 * @param xy Used to set the gene expressions x and protein concentrations y (if translation)
	 * @param dxydt Variations of x and y (output)
	 */
	public void computeDxydt(double[] xy, double[] dxydt) {
		
		boolean modelTranslation = GnwSettings.getInstance().getModelTranslation();
		//double[] dxydt = new double[xy.length];
		int size = getSize();
		
		for (int i=0; i<size; i++)
			x_.set(i, xy[i]);
		
		if (modelTranslation)
			for (int i=0; i<size; i++)
				y_.set(i, xy[size+i]);
		else
			y_ = x_;
		
		// dxydt temporarily used to store the production rates of mRNA
		computeMRnaProductionRates(dxydt);
		
		for (int i=0; i < size; i++)
		{
			if ( nodes_.get(i) instanceof Gene)
			{
				dxydt[i] = dxydt[i] - ((Gene)nodes_.get(i)).computeMRnaDegradationRate(x_.get(i));
			}
			else if (nodes_.get(i) instanceof Environment)
			{
				dxydt[i] = ((Environment)nodes_.get(i)).computeDxdt(time_);
			}		
		}
		
		if (modelTranslation)
		{
			for (int i=0; i<size; i++)
			{
				if ( nodes_.get(i) instanceof Gene)
				{
					dxydt[size+i] = ((Gene)nodes_.get(i)).getMaxTranslation()*x_.get(i) - ((Gene)nodes_.get(i)).computeProteinDegradationRate(y_.get(i));
				}
				else if (nodes_.get(i) instanceof Environment)
				{
					dxydt[size+i] = ((Environment)nodes_.get(i)).computeDxdt(time_);
				}	
			}
		}
	}

	
	// ----------------------------------------------------------------------------
	
	/** 
	 * Load a gene network from an SBML file. Overrides Structure.load(). Format must
	 * be equal GeneNetwork.SBML. Note, the SBML file must be in the exact same format
	 * as the SBML files produced by writeSBML(). In particular, we assume that reactions are listed
	 * *ordered* as we do in writeSBML().
	 * @param filename URL to the file describing the network to load
	 * @param format File format (GML, DOT, etc.)
	 * @throws IOException 
	 */
	@SuppressWarnings("unchecked")
	public void load(URL filename, int format) throws IOException {
		
		if (format != SBML)
			throw new IllegalArgumentException("GeneNetwork.load(): format must be equal GeneNetwork.SBML");
		
		// ----------------------
		// Load the SBML document
		SBMLLevel2Document sbml2Doc = null;
		Model model = null;
		
		File f = new File(filename.getPath());
		id_ = f.getName();
		
		URLConnection uc = filename.openConnection();

		sbml2Doc = SBMLLevel2Document.readDocument(uc.getInputStream());
		model = sbml2Doc.getModel();
		
		List<Species> species = model.getSpecies();
		int size = species.size();

		// -----------------------------------------
		// Set the network size and create the genes
		// do not count the species _void_ and _time_
		setSize(species.size()-2);

		for (int g=0; g < size; g++) 
		{
			String id = species.get(g).getId();
			Node node = GrnNodeFactory.createNode(id, this);
			if( node != null)
			{
				node.setLabel(species.get(g).getId());
				nodes_.add(node);
			}
		}
		
		ArrayList<String> parameterNames = new ArrayList<String>(); // the names of the parameters
		ArrayList<Double> parameterValues = new ArrayList<Double>(); // the values of the parameters
		ArrayList<GrnNode> inputIndexes = new ArrayList<GrnNode>(); // the indexes of the inputs
		
		// 2 loops for one gene: both synthesis and degradation reactions
		// (we assume that reactions are listed *ordered* as we do in writeSBML())
		for (
		Iterator<?> i = model.getReactions().iterator(); i.hasNext(); ) {
			
			// get the next reaction and the modifiers
			Reaction r = (Reaction) i.next();
			GrnNode target = (GrnNode)getNode(getGeneReactantId(r));
			
			if (target == null && !getGeneReactantId(r).equals("_time_"))
				log_.log(Level.WARNING, "Error: " + getGeneReactantId(r));
	
			// get all the genes modifying the target gene
			List<?> modifiers = r.getModifier();

			// Set edges (if this is a degradation reaction, there are no modifiers and this loop is not entered)
			for (
			Iterator<?> j = modifiers.iterator(); j.hasNext();) {
				ModifierSpeciesReference speciesReference = (ModifierSpeciesReference) j.next();
				GrnNode source = (GrnNode)getNode(speciesReference.getSpecies());
				if(source != null )
				{
					inputIndexes.add(source);
					// The edge type is unknown for now, it is initialized later
					edges_.add(new Edge(source, target, Edge.UNKNOWN));
				}
			}
		
			// compile parameters list of gene
			for (
			Iterator<?> k = r.getKineticLaw().getParameter().iterator(); k.hasNext();) {
				Parameter param = (Parameter)k.next();
				parameterNames.add(param.getId());
				parameterValues.add(param.getValue());
			}
			
			// in the second iteration for this gene
			if (r.getId().contains("degradation") || target instanceof Environment)
			{
				// set parameters in gene
				GrnNode g = (GrnNode)target;
				g.initialization(parameterNames, parameterValues, inputIndexes);
				inputIndexes = new ArrayList<GrnNode>(); // don't clear because the reference was copied to the gene
				parameterNames.clear(); // reset (they were not copied)
				parameterValues.clear();
			}
		}
		setEdgeTypesAccordingToDynamicalModel();
		signed_ = true;
	}


	// ----------------------------------------------------------------------------

	/**
	 * Set the type (enhancing, inhibiting, dual, or unknown) of all edges according
	 * to the dynamically model.
	 */
	public void setEdgeTypesAccordingToDynamicalModel() {
		
		for (int i=0; i<nodes_.size(); i++) 
		{
			if (nodes_.get(i) instanceof Gene)
			{
				((Gene)nodes_.get(i)).setInputEdgeTypesAccordingToDynamicalModel();			
			}
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/**
	 * Return a list of all genes that perturbe the gene given in parameter.
	 * @param gene Gene perturbed
	 * @return Instance of all genes that perturbed the gene given in parameter.
	 */
	public ArrayList<GrnNode> getInputs(Node node) {
		
		ArrayList<GrnNode> inputs = new ArrayList<GrnNode>();

		for (int i=0; i < edges_.size(); i++) {
			if (edges_.get(i).getTarget().getLabel().equals(node.getLabel())) {
				String source = edges_.get(i).getSource().getLabel();
				inputs.add((GrnNode) getNode(source));
			}
		}
		
		return inputs;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Return gene id of a reaction.
	 * @param r Instance of Reaction
	 * @return Gene id of the reaction r
	 */
	public String getGeneReactantId(Reaction r) {
		
		// Example: G1_synthesis or G1_degradation -> id = "G1"
		//StringTokenizer st = new StringTokenizer(r.getId(), "_");
		//String id = st.nextToken();
		// problem: if the id also has '_', the above code fails
		String id = r.getId();
		id = id.substring(0, id.lastIndexOf("_"));
		return id;
	}
	
	
	// ----------------------------------------------------------------------------

	/**
	 * Save the gene network to an SBML file. If the argument is null, use the network id.
	 * @param filename URL to the file describing the network to load
	 * @throws IOException
	 */
	@SuppressWarnings({ "unused"})
	public int writeSBML(URL filename) throws IOException {
		
		Model model = new Model(id_);
		model.getNotes ().add (comment_); // save network description
		
		log_.log(Level.INFO, "Writing file " + filename.getPath());
		
		int size = getSize();
		
		Compartment compartment = new Compartment ("cell", "cell");
		compartment.setSize (1);
		model.addCompartment (compartment);

		// Hack in pattern expression in sbml2 library to allow any char to compose
		// the identifier of a gene (hack in top of SBaseId.java)
		Species[] species = new Species[size+2];
  
		for (int s=0; s < size; s++) { // save gene as species
//			species[s] = new Species(nodeIds_.get(s), nodeIds_.get(s));
			species[s] = new Species(nodes_.get(s).getLabel(), nodes_.get(s).getLabel());
			species[s].setCompartment(compartment);
			//species[s].setInitialAmount(?); // maybe save the wild-type steady state?
			model.addSpecies(species[s]);
		}
		
		// create the void species
		species[size] = new Species("_void_", "_void_");
		species[size].setCompartment(compartment);
		species[size].setInitialAmount(0);
		species[size].setBoundaryCondition(true);
		species[size].setConstant(true);
		model.addSpecies(species[size]);
		
		// create the time species
		species[size+1] = new Species("_time_", "_time_");
		species[size+1].setCompartment(compartment);
		species[size+1].setInitialAmount(0);
		model.addSpecies(species[size+1]);
				
		// SET SYNTHESIS AND DEGRADATION REACTIONS FOR EVERY GENE
		boolean environment = false;
		for (int i=0; i<size; i++) {

			// the ID of gene i
//			String currentGeneID = nodeIds_.get(i);
			Node current = nodes_.get(i);
			String currentGeneID = current.getLabel();
			// The modifiers (regulators) of gene i
			if ( current instanceof Gene)
			{			
				ArrayList<GrnNode> inputIndexes = ((GrnNode)current).getInputGenes();
  
				// SYNTHESIS REACTION
				KineticLaw kineticLaw = new KineticLaw ();
				String reactionId = currentGeneID + "_synthesis";
				String name = reactionId + ": " + ((GrnNode)current).toString();
				Reaction reaction = new Reaction (reactionId, name);
				reaction.setReversible (false);
				reaction.addReactant(species[size]);
				reaction.addProduct(species[i]);
				for (int r=0; r<inputIndexes.size(); r++) // set gene modifiers
//					reaction.addModifier(species[inputIndexes.get(r)]);
					reaction.addModifier(species[nodes_.indexOf(inputIndexes.get(r))]);
  
				ArrayList<String> names = new ArrayList<String>(); // parameters names
				ArrayList<Double> values = new ArrayList<Double>(); // parameters values
				((GrnNode)nodes_.get(i)).compileParameters(names, values);
			
				// save gene parameters (note, the first param is the degradation rate)
				for (int p=1; p<names.size(); p++) {
					String id = names.get(p);// + "_" + currentGeneID;
					kineticLaw.addParameter(new Parameter(id, id, values.get(p)));
				}
				String mathSynth = createSynthesisFormula((HillGene)current);
				kineticLaw.setMath(mathSynth);

				reaction.setKineticLaw(kineticLaw);
				model.addReaction(reaction);
  
				// DEGRADATION REACTION
				kineticLaw = new KineticLaw ();
				reactionId = currentGeneID + "_degradation";
				reaction = new Reaction (reactionId, reactionId);
				reaction.setReversible(false);
				reaction.addReactant(species[i]);
				reaction.addProduct(species[size]);
  
				String id = names.get(0);
				kineticLaw.addParameter(new Parameter(id, id, values.get(0)));
				String mathDeg = createDegradationFormula(current);
				kineticLaw.setMath(mathDeg);
				reaction.setKineticLaw (kineticLaw);
				model.addReaction (reaction);
			}
			else if(current instanceof Environment)
			{
				if (!environment)
				{
					//We don't have an environment yet Create a "time" node to ensure we can use time explicitly.
					Environment dummy = new LinearEnvironment(this, new int[] {}, "_time_", 0, 1);
					ArrayList<String> names = new ArrayList<String>(); // parameters names
					ArrayList<Double> values = new ArrayList<Double>(); // parameters values
					dummy.compileParameters(names, values);
					dummy.writeSBML(model, species[size+1], species[size], names, values);
					environment = true;
				}
				ArrayList<String> names = new ArrayList<String>(); // parameters names
				ArrayList<Double> values = new ArrayList<Double>(); // parameters values
				((GrnNode)nodes_.get(i)).compileParameters(names, values);

				((Environment) current).writeSBML(model, species[i], species[size], names, values);
			}
		}
		
		try
		{
			// PRINT FILE
			SBMLLevel2Document sbml2Doc = new SBMLLevel2Document(model);
			if (filename == null)
				filename = new URL(id_ + ".xml");
			sbml2Doc.writeDocument(new FileWriter(new File(filename.toURI()), false));
		}
		catch (URISyntaxException e)
		{
			log_.log(Level.WARNING, "URISyntaxException", e);
		}
		
		return 0;
	}
	
	/**
	 * Create the MATHML strings for all regulatory modules
	 * @param gene the gene to create the MathML for
 	 * @return Per Module, an ArrayList<String> that represents all lines of the MathML object
	 */
	private ArrayList<ArrayList<String>> createModulesFormula(HillGene gene)
	{
		//We create his formula as a vector of strings, such that we can fix the tabs later on
		ArrayList<ArrayList<String>> moduleFormulas = new ArrayList<ArrayList<String>>();
		int nextInput = 0;

		for (int i = 0; i < gene.getRegulatoryModules().size(); i++)
		{
			int currInput = nextInput;
			RegulatoryModule module = gene.getRegulatoryModules().get(i);
			ArrayList<String> moduleFormula = new ArrayList<String>();
			// get the inputs for this module
			moduleFormula.add("<apply><divide/>");

			//get the formula for the numerator
			if (module.getNumActivators() == 0)
			{
				moduleFormula.add("\t<mn>1</mn>");
			}
			else if (module.getNumActivators() == 1)
			{
				String geneName = gene.getInputGenes().get(currInput++).getLabel();
				createActivationFormula(currInput, moduleFormula, geneName, 0);
			}
			else
			{
				moduleFormula.add("\t<apply><times/>");

				for (int j=0; j<module.getNumActivators(); j++)
				{
					String geneName = gene.getInputGenes().get(currInput++).getLabel();
					createActivationFormula(currInput, moduleFormula, geneName, 1);
				}
				moduleFormula.add("\t</apply>");
			}
			currInput = nextInput;
			//get the formula for the denominator

			if( module.getNumInputs() == 1)
			{
				String geneName = gene.getInputGenes().get(currInput++).getLabel();
				moduleFormula.add("\t<apply><plus/>");
				moduleFormula.add("\t\t<cn>1</cn>");
				createActivationFormula(currInput, moduleFormula, geneName, 1);
				moduleFormula.add("\t</apply>");
			}
			else
			{
				moduleFormula.add("\t<apply><plus/>");
				moduleFormula.add("\t\t<apply><times/>");
				for (int j=0; j<module.getNumInputs(); j++)
				{
					String geneName = gene.getInputGenes().get(currInput++).getLabel();
					moduleFormula.add("\t\t\t<apply><plus/>");
					moduleFormula.add("\t\t\t\t<cn>1</cn>");
					createActivationFormula(currInput, moduleFormula, geneName, 3);
					moduleFormula.add("\t\t\t</apply>");
				}
				moduleFormula.add("\t\t</apply>");
				moduleFormula.add("\t</apply>");
			}
			moduleFormula.add("</apply>");
			 
			moduleFormulas.add(moduleFormula);
			nextInput = currInput;
		}
		
		return moduleFormulas;
	}

	/**
	 * Append the MATHML code to the current object
	 * @param geneIndx index of the gene
	 * @param moduleFormula MathML object to append to
	 * @param geneName name of the gene
	 */
	private void createActivationFormula(int geneIndx, ArrayList<String> moduleFormula, String geneName, Integer numLeadingTabs)
	{
		//appends the mathML for the activation factor of this gene
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t<apply><power/>");
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t\t<apply><divide/>");
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t\t\t<ci>" + geneName + "</ci>");
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t\t\t<ci>k_" + geneIndx + "</ci>");
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t\t</apply>");
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t\t<ci>n_" + geneIndx + "</ci>");
		moduleFormula.add("\t".repeat(numLeadingTabs) + "\t</apply>");
	}
	
	/**
	 * Create the MATHML string for the synthesis reaction of a gene
	 * @param gene the gene in question
 	 * @return The MathML string
	 */
	private String createSynthesisFormula(HillGene gene) 
	{
		//This has to follow the logic of HillgGene.computeMRnaProductionRate precisely
		String math = "\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
		ArrayList<ArrayList<String>> mathModules = createModulesFormula(gene);
		
		if( gene.getRegulatoryModules().size() == 0)
		{
			math += "\t\t\t<apply><times/>\n";
			math += "\t\t\t\t<ci> " + "cell" + "</ci>\n";
			math += "\t\t\t\t<ci>max</ci>\n";
			math += "\t\t\t\t<ci>a_0</ci>\n";
			math += "\t\t\t</apply>\n";
			math += "\t\t</math>\n";
		}
		else
		{
			math += "\t\t\t<apply><times/>\n";
			math += "\t\t\t\t<ci> " + "cell" + "</ci>\n";
			math += "\t\t\t\t<ci>max</ci>\n";
			math += "\t\t\t\t<apply><plus/>\n";
		
			for (int i=0; i<Math.pow(2,gene.getRegulatoryModules().size()); i++) 
			{
				// binary representation of state i
				String s = Integer.toBinaryString(i); // note, leading 0's are not included
			
				math += "\t\t\t\t\t<apply><times/>\n";
				math += "\t\t\t\t\t\t<ci>a_" + i + "</ci>\n";

				// if module j is active in this state, multiply with m_j, otherwise with (1-m_j)
				for (int j=0; j<gene.getRegulatoryModules().size(); j++) 
				{
					// if s.length()-j-1 smaller than zero, the corresponding module is off (it's one of the leading zeros)
					if (s.length()-j-1 >= 0 && s.charAt(s.length()-j-1) == '1')
					{
						for (String line : mathModules.get(j))
						{
							math += "\t\t\t\t\t\t" + line + "\n";
						}
					}
					else
					{
						math += "\t\t\t\t\t\t<apply><minus/>\n";
						math += "\t\t\t\t\t\t\t<cn>" + 1 + "</cn>\n";
						for (String line : mathModules.get(j))
						{
							math += "\t\t\t\t\t\t\t" + line + "\n";
						}
						math += "\t\t\t\t\t\t</apply>\n";
					}
				}
				math += "\t\t\t\t\t</apply>\n";
			}
			math += "\t\t\t\t</apply>\n";
			math += "\t\t\t</apply>\n";

			math += "\t\t</math>\n";
		}

		return math;
	}

	/**
	 * Create the MATHML string for the degradation reaction of a gene
	 * @param node the node of the gene in question
 	 * @return The MathML string
	 */
	private String createDegradationFormula(Node node) 
	{
		String math = "\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n"
				+ "\t\t\t<apply>\n"
				+ "\t\t\t\t<times/>\n"
				+ "\t\t\t\t\t<ci> " + "cell" + "</ci>\n"
				+ "\t\t\t\t\t<ci> " + "delta" + "</ci>\n"
				+ "\t\t\t\t\t<ci> " + node.getLabel() + "</ci>\n"
				+ "\t\t\t</apply>\n"
				+ "\t\t</math>\n";
		return math;
	}
	
	// ----------------------------------------------------------------------------
	
	/** Get a string with all gene labels separated by tabs (often used as header when printing data) */
	public String getHeader(boolean includeProteins) {

		String header = nodes_.get(0).getLabel();
		for (int i=1; i<nodes_.size(); i++)
			header += "\t" + nodes_.get(i).getLabel();
		
		if (includeProteins)
			for (int i=1; i<nodes_.size(); i++)
				header += "\t" + nodes_.get(i).getLabel() + "_prot";
		header += "\n";
		return header;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get the header as an array of String */
	public String[] getHeaderArray() {

		String[] header = new String[nodes_.size()];
		
		for (int i=0; i<nodes_.size(); i++)
			header[i] = nodes_.get(i).getLabel();
		
		return header;
	}

	
	// ----------------------------------------------------------------------------
	
	/**
	 * Create a subnetwork from a list of nodes.
	 * @param name Name of the new sub-network
	 * @param seed Group of nodes
	 * @return Sub-network
	 */
	public GeneNetwork getSubnetwork(String name, ArrayList<Node> seed)
	{
		GeneNetwork network = new GeneNetwork();
		network.id_ = name;
	
		// Copies of all the nodes
		for (int i=0; i < seed.size(); i++) {
			Node gene = seed.get(i);
			if (gene instanceof HillGene) {
//				System.out.println("GeneHill");
				network.addNode( ((HillGene) gene).copy() );
			}
		}
		
		// Set the topology (edges)
		Node source, target;
		Edge edge;
		
		for (int i=0; i < edges_.size(); i++)
		{
			edge = edges_.get(i);
			// Search in network if the source and target nodes exist
			source =  network.getNode( edge.getSource().getLabel() );
			target =  network.getNode( edge.getTarget().getLabel() );
			
			if (source != null && target != null)
				network.addEdge( new Edge(source, target, edge.getType()) ); //, edge.getLabel()) );
		}
		return network;
	}	
		
	
	// ----------------------------------------------------------------------------

	/** Set the flags Gene.isTf_ of all regulators, set numTfs_ */
	public void markTfs() {
		
		for (int i=0; i<edges_.size(); i++)
			((Gene) edges_.get(i).getSource()).setIsTf(true);

		numTfs_ = 0;
		for (int i=0; i<nodes_.size(); i++)
			if (isGene(i) && getGene(i).getIsTf())
				numTfs_++;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the index of a randomly selected TF */
	public int getRandomTf() {
		
		if (numTfs_ < 1)
			throw new RuntimeException("There are no TFs in the network (make sure markTfs() has been called)");
		
		int tfIndex = GnwSettings.getInstance().getUniformDistribution().nextIntFromTo(0, numTfs_-1);
		int tfCount = -1;
		int index = 0;

		for (; index<nodes_.size(); index++) 
		{
			if (isGene(index) && getGene(index).getIsTf()) {
				tfCount++;
				if (tfCount == tfIndex)
					break;
			}
		}
		
		assert index < nodes_.size();
		return index;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Add the give fraction of random edges to the network (outgoing edges are only added to existing regulators, no duplicate edges will be created) */
	public void addRandomEdges(double fractionRandomEdges) {
		
		int numRandomEdges = (int) (fractionRandomEdges * edges_.size());
		Uniform uni = GnwSettings.getInstance().getUniformDistribution();
		boolean noSelfLoops = GnwSettings.getInstance().getIgnoreAutoregulatoryInteractionsInEvaluation();
		
		int counter = 0;
		while (counter < numRandomEdges) {
			
			Node tf = nodes_.get(getRandomTf());
			Node gene = nodes_.get( uni.nextIntFromTo(0, nodes_.size()-1) );
			
			if (noSelfLoops && tf == gene)
				continue;
			
			boolean duplicateEdge = false;
			for (int j=0; j<edges_.size(); j++) {
				if (edges_.get(j).getSource() == tf && edges_.get(j).getTarget() == gene) {
					duplicateEdge = true;
					break;
				}
			}
			if (duplicateEdge)
				continue;
			
			addEdge(new Edge(tf, gene));
			counter++;
		}
	}

	
	// ============================================================================
	// PRIVATE FUNCTIONS
		
	/**
	 * Compute the production rate for all the genes. This function is identical to
	 * the one below, except that the argument is passed as double[] instead of
	 * DoubleMatrix1d. We need both because the SDE implementation works with
	 * DoubleMatrix1d and the ODE implementation works with arrays.
	 * NOTE: x_ must be set before calling this function.
	 * @param productionRates Return the production rates of all the genes
	 */
	protected void computeMRnaProductionRates(double[] productionRates) {
		
		int size = getSize();
		
		for (int i=0; i < size; i++) 
		{
			if ( nodes_.get(i) instanceof Gene)
			{
				if (GnwSettings.getInstance().getModelTranslation())
					productionRates[i] = ((Gene)nodes_.get(i)).computeMRnaProductionRate(i, y_);
				else
					productionRates[i] = ((Gene)nodes_.get(i)).computeMRnaProductionRate(i, x_);
			}
		}
	}
	

	/**
	 * Compute the production rate for all the genes. This function is identical to
	 * the one above, except that the argument is passed as DoubleMatrix1d instead of
	 * double[]. We need both because the SDE implementation works with
	 * DoubleMatrix1d and the ODE implementation works with arrays.
	 * NOTE: x_ must be set before calling this function.
	 * @param productionRates Return the production rates of all the genes
	 */
	protected void computeMRnaProductionRates(DoubleMatrix1D productionRates) {
		
		int size = getSize();
		
		for (int i=0; i < size; i++) {
			if (GnwSettings.getInstance().getModelTranslation())
			{
				if(isGene(i))
				{
					productionRates.set(i, ((Gene)nodes_.get(i)).computeMRnaProductionRate(i, y_));
				}
				else
				{
					productionRates.set(i, ((Environment)nodes_.get(i)).computeDxdt(time_));
				}
			}
			else
			{
				if(isGene(i))
				{
					productionRates.set(i, ((Gene)nodes_.get(i)).computeMRnaProductionRate(i, x_));
				}
				else
				{
					productionRates.set(i, ((Environment)nodes_.get(i)).computeDxdt(time_));
				}
			}
		}
	}
	
	
	// ----------------------------------------------------------------------------
	
	protected void computeProteinProductionRates(DoubleMatrix1D productionRates) {
		
		int size = getSize();
		
		for (int i=0; i<size; i++)
		{
			if(isGene(i))
			{
				productionRates.set(i, ((Gene)nodes_.get(i)).computeProteinProductionRate(x_.get(i)));			
			}
			else
			{
				productionRates.set(i, ((Environment)nodes_.get(i)).computeDxdt(time_));
			}
		}
	}
	
	
	// ----------------------------------------------------------------------------
	
	protected void computeMRnaDegradationRates(DoubleMatrix1D degradationRates) {
		
		int size = getSize();
		
		for (int i=0; i<size; i++)
		{
			if(isGene(i))
			{
				degradationRates.set(i, ((Gene)nodes_.get(i)).computeMRnaDegradationRate(x_.get(i)));
			}
			else
			{
				degradationRates.set(i,0);
			}
		}
	}
	
	
	// ----------------------------------------------------------------------------
	
	protected void computeProteinDegradationRates(DoubleMatrix1D degradationRates) {
		
		int size = getSize();
		
		for (int i=0; i<size; i++)
		{
			if(isGene(i))
			{
				degradationRates.set(i, ((Gene)nodes_.get(i)).computeProteinDegradationRate(y_.get(i)));
			}
			else
			{
				degradationRates.set(i,0);
			}
		}
	}
	
	public void addEnvironment(Environment env)
	{
		nodes_.add(env);
		for (int indx : env.getTargetGenes())
		{
			Node target = nodes_.get(indx);
			edges_.add( new Edge(env, target, Edge.ENHANCER) );
		}
	}
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public Gene getGene(int i) { return (Gene)nodes_.get(i); }
	public boolean isGene(int i) { return (nodes_.get(i) instanceof Gene);}
	public DoubleMatrix1D getX() { return x_; }
	public void setX(DoubleMatrix1D x) { x_ = x; }
	
	public DoubleMatrix1D getY() { return y_; }
	public void setY(DoubleMatrix1D y) { y_ = y; }
	
	public int getNumTfs() { return numTfs_; }
	
	public void increaseTime(double delta) { time_ += delta;}
	
	public void resetTime() { time_ = 0;}
	
	public double getTime() { return time_;}
}
