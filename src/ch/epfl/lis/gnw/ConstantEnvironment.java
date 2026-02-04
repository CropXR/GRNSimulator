package ch.epfl.lis.gnw;

import java.util.ArrayList;

import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.Species;

public class ConstantEnvironment extends Environment
{
	private double envValue_ = 0;
	
	/** Default constructor */
	public ConstantEnvironment() 
	{
		super();
	}
	
	public ConstantEnvironment(GeneNetwork grn)
	{
		super(grn);
	}
	
	public ConstantEnvironment(GeneNetwork grn, int[] targetGenes, String label, double value) 
	{
		super(grn, targetGenes, label);
		envValue_ = value;
	}

	@Override
	public double computeActivity()
	{
		return envValue_;
	}

	@Override
	public void compileParameters(ArrayList<String> names, ArrayList<Double> values) 
	{
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setInputEdgeTypesAccordingToDynamicalModel()
	{
		return;	
	}

	@Override
	public void writeSBML(Model model, Species thisSpecies, Species voidSpecies) 
	{
		//Constant environment, so only initial condition needed
		thisSpecies.setInitialAmount(envValue_);
		thisSpecies.setConstant(true);
		
		//add math string because it apparently required
		String math = "\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
		math += "\t\t\t<cn>0</cn>\n";
		math += "\t\t</math>\n";

		KineticLaw kineticLaw = new KineticLaw ();
		String reactionId = label_ + "_synthesis";
		Reaction reaction = new Reaction (reactionId, reactionId);
		reaction.setReversible(false);
		reaction.addReactant(voidSpecies);
		reaction.addProduct(thisSpecies);
				
		kineticLaw.addParameter(new Parameter("envValue", "envValue", envValue_));

		kineticLaw.setMath(math);

		reaction.setKineticLaw(kineticLaw);
		model.addReaction(reaction);
	
		return;	
	}

	@Override
	public double computeDxdt() 
	{
		return 0;
	}

	@Override
	public double getInitialCondition() 
	{
		return envValue_;
	}

	@Override
	public void initialization(ArrayList<String> paramNames, ArrayList<Double> paramValues,
			ArrayList<GrnNode> inputGenes) 
	{
		envValue_ = paramValues.get(paramNames.indexOf("envValue"));
		
		return;
	}	
}
