package ch.epfl.lis.gnw;

import java.util.ArrayList;

import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.Species;

public class LinearEnvironment extends Environment
{
	private double startVal_ = 0;
	private double rate_ = 0;

	/** Default constructor */
	public LinearEnvironment() 
	{
		super();
	}
	
	public LinearEnvironment(GeneNetwork grn)
	{
		super(grn);
	}
	
	public LinearEnvironment(GeneNetwork grn, int[] targetGenes, String label, double start, double rate) 
	{
		super(grn, targetGenes, label);
		startVal_ = start;
		rate_ = rate;
	}

	@Override
	public double computeActivity(double time)
	{
		return startVal_ + time * rate_;
	}

	@Override
	public void compileParameters(ArrayList<String> names, ArrayList<Double> values) 
	{
		names.add("startVal");
		values.add(startVal_);	
		
		
		names.add("rate");
		values.add(rate_);		
	}

	@Override
	public void writeSBML(Model model, Species thisSpecies, Species voidSpecies, ArrayList<String> names, ArrayList<Double> values) 
	{
		thisSpecies.setInitialAmount(startVal_);
		
		KineticLaw kineticLaw = new KineticLaw ();
		String math = "\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
		math += "\t\t\t<ci>rate</ci>\n";
		math += "\t\t</math>\n";

		String reactionId = label_ + "_synthesis";
		Reaction reaction = new Reaction (reactionId, reactionId);
		reaction.setReversible(false);
		reaction.addReactant(voidSpecies);
		reaction.addProduct(thisSpecies);
		
		kineticLaw.setMath(math);
		
		for (int p=0; p<names.size(); p++) 
		{
			String id = names.get(p);// + "_" + currentGeneID;
			kineticLaw.addParameter(new Parameter(id, id, values.get(p)));
		}

		reaction.setKineticLaw(kineticLaw);
		model.addReaction(reaction);
	
		return;	
	}

	@Override
	public double computeDxdt(double time) 
	{
		return rate_;
	}

	@Override
	public double getInitialCondition()
	{
		return startVal_;
	}

	@Override
	public void initialization(ArrayList<String> paramNames, ArrayList<Double> paramValues,
			ArrayList<GrnNode> inputGenes) 
	{
		rate_ = paramValues.get(paramNames.indexOf("rate"));
		startVal_ = paramValues.get(paramNames.indexOf("startVal"));
		
		return;
	}	
}