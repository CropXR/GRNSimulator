package ch.epfl.lis.gnw;

import java.util.logging.Logger;

import jigcell.sbml2.Model;
import jigcell.sbml2.Species;

abstract public class Environment extends GrnNode
{	
	/** Reference to the gene regulatory network of this environment */
	protected double time_ = 0;
	private int[] targetGenes_; 
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor */
	public Environment() 
	{
		super();
	}
	
	public Environment(GeneNetwork grn)
	{
		super();
		grn_ = grn;
		targetGenes_ = null;
		label_ = "";
	}
	
	public Environment(GeneNetwork grn, int[] targetGenes, String label) 
	{
		super();
		grn_ = grn;
		targetGenes_ = targetGenes;
		label_ = label;
	}	
	
	abstract public double computeActivity();
	
	public void increaseTime(double delta)
	{
		time_ += delta;
	}
	
	public int[] getTargetGenes()
	{
		return targetGenes_;
	}
	
	abstract public void writeSBML(Model model, Species thisSpecies, Species voidSpecies);
	abstract public double computeDxdt();
	abstract public double getInitialCondition();

	/** Logger for this class */
	@SuppressWarnings("unused")
	private static Logger log_ = Logger.getLogger(Environment.class.getName());
}
