package ch.epfl.lis.gnw;

import java.util.ArrayList;
import java.util.logging.Logger;

import ch.epfl.lis.networks.Node;

abstract public class GrnNode extends Node
{
    /** Logger for this class */
	@SuppressWarnings("unused")
	private static Logger log_ = Logger.getLogger(GrnNode.class.getName());
	
	/** Reference to the gene regulatory network of this gene */
	protected GeneNetwork grn_ = null;
	
	/** 
	 * An array with references to the inputs of this gene. The order defines which
	 * input goes where. E.g., inputGenes_[0] is the first input of the first module. 
	 */
	protected ArrayList<GrnNode> inputGenes_ = null;
	
	/** Default constructor */
	public GrnNode() {
		super();
	}
	
	
	/** Constructor, sets the given gene network */
	public GrnNode(GeneNetwork grn) {
		super();
		grn_ = grn;
	}	

	public ArrayList<GrnNode> getInputGenes() {return inputGenes_; }
	public void setInputGenes(ArrayList<GrnNode> inputs) { inputGenes_ = inputs; }
	
	abstract public void compileParameters(ArrayList<String> names, ArrayList<Double> values);
	abstract public void initialization(ArrayList<String> paramNames, ArrayList<Double> paramValues, ArrayList<GrnNode> inputGenes);

}
