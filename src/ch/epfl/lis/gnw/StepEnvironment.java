package ch.epfl.lis.gnw;

import java.util.ArrayList;

import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.Species;

public class StepEnvironment extends Environment {

	private double firstEnvValue_ = 0;
	private double secondEnvValue_ = 0;
	private double stepTime_ = 0;
	private double steepness_ = 0;
	/** Default constructor */
	public StepEnvironment() 
	{
		super();
	}
	
	public StepEnvironment(GeneNetwork grn)
	{
		super(grn);
	}
	
	public StepEnvironment(GeneNetwork grn, int[] targetGenes, String label, 
			double firstEnvValue, double secondEnvValue, double stepTime, double steepness ) 
	{
		super(grn, targetGenes, label);
		firstEnvValue_ = firstEnvValue;
		secondEnvValue_ = secondEnvValue;
		stepTime_ = stepTime;
		steepness_ = steepness;
	}
	
	@Override
	public double computeActivity(double time) 
	{
		return firstEnvValue_ + (secondEnvValue_ - firstEnvValue_) * Math.pow(time, steepness_) /
				(Math.pow(stepTime_, steepness_) + Math.pow(time, steepness_));
		
	}

	@Override
	public void writeSBML(Model model, Species thisSpecies, Species voidSpecies, ArrayList<String> names,
			ArrayList<Double> values) 
	{
		thisSpecies.setInitialAmount(firstEnvValue_);
		
		KineticLaw kineticLaw = new KineticLaw ();
		String math = "\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
		math += "\t\t<apply><divide/>\n";
		math += "\t\t\t<apply><times/>\n";
		math += "\t\t\t\t<ci>steepness</ci>\n";
		math += "\t\t\t\t<apply><power/>\n";
		math += "\t\t\t\t\t<ci>_time_</ci>\n";
		math += "\t\t\t\t\t<apply><minus/>\n";
		math += "\t\t\t\t\t\t<ci>steepness</ci>\n";
		math += "\t\t\t\t\t\t<cn>1</cn>\n";
		math += "\t\t\t\t\t</apply>\n";
		math += "\t\t\t\t</apply>\n";
		math += "\t\t\t\t<apply><power/>\n";
		math += "\t\t\t\t\t<ci>stepTime</ci>\n";
		math += "\t\t\t\t\t<ci>steepness</ci>\n";
		math += "\t\t\t\t</apply>\n";
		math += "\t\t\t\t<apply><minus/>\n";
		math += "\t\t\t\t\t<ci>secondEnvValue</ci>\n";
		math += "\t\t\t\t\t<ci>firstEnvValue</ci>\n";
		math += "\t\t\t\t</apply>\n";
		math += "\t\t\t</apply>\n";
		math += "\t\t\t<apply><power/>\n";
		math += "\t\t\t\t<apply><plus/>\n";
		math += "\t\t\t\t\t<apply><power/>\n";
		math += "\t\t\t\t\t\t<ci>_time_</ci>\n";
		math += "\t\t\t\t\t\t<ci>steepness</ci>\n";
		math += "\t\t\t\t\t</apply>\n";
		math += "\t\t\t\t\t<apply><power/>\n";
		math += "\t\t\t\t\t\t<ci>stepTime</ci>\n";
		math += "\t\t\t\t\t\t<ci>steepness</ci>\n";
		math += "\t\t\t\t\t</apply>\n";
		math += "\t\t\t\t</apply>\n";
		math += "\t\t\t\t<cn>2</cn>";
		math += "\t\t\t</apply>\n";
		math += "\t\t</apply>\n";
		math += "\t\t</math>\n";

		String reactionId = label_ + "_synthesis";
		Reaction reaction = new Reaction (reactionId, reactionId);
		reaction.setReversible(false);
		reaction.addReactant(voidSpecies);
		reaction.addProduct(thisSpecies);
		reaction.addModifier(new Species("_time_", "_time_"));

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
		double denominatorpart = Math.pow(time, steepness_) + Math.pow(stepTime_, steepness_);
		return steepness_ * Math.pow(time, steepness_ - 1) * Math.pow(stepTime_, steepness_) * (secondEnvValue_ - firstEnvValue_) 
				/ (denominatorpart * denominatorpart);
	}

	@Override
	public double getInitialCondition() 
	{
		return firstEnvValue_;
	}

	@Override
	public void compileParameters(ArrayList<String> names, ArrayList<Double> values) 
	{
		names.add("firstEnvValue");
		values.add(firstEnvValue_);		
		names.add("secondEnvValue");
		values.add(secondEnvValue_);
		names.add("stepTime");
		values.add(stepTime_);
		names.add("steepness");
		values.add(steepness_);
	}

	@Override
	public void initialization(ArrayList<String> paramNames, ArrayList<Double> paramValues,
			ArrayList<GrnNode> inputGenes) 
	{
		firstEnvValue_ = paramValues.get(paramNames.indexOf("firstEnvValue"));
		secondEnvValue_ = paramValues.get(paramNames.indexOf("secondEnvValue"));
		stepTime_ = paramValues.get(paramNames.indexOf("stepTime"));
		steepness_ = paramValues.get(paramNames.indexOf("steepness"));
		return;
	}
}
