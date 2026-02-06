package ch.epfl.lis.gnw;

import java.util.ArrayList;

import jigcell.sbml2.KineticLaw;
import jigcell.sbml2.Model;
import jigcell.sbml2.Parameter;
import jigcell.sbml2.Reaction;
import jigcell.sbml2.Species;

public class OscillatingEnvironment extends Environment
{
	private double phase_ = 0;
	private double frequency_ = 0;
	private double amplitude_ = 0;
	private double offset_ = 0;

	/** Default constructor */
	public OscillatingEnvironment() 
	{
		super();
	}
	
	public OscillatingEnvironment(GeneNetwork grn)
	{
		super(grn);
	}
	
	public OscillatingEnvironment(GeneNetwork grn, int[] targetGenes, String label, double phase, double frequency, double amplitude, double offset) 
	{
		super(grn, targetGenes, label);
		phase_ = phase;
		frequency_ = frequency;
		amplitude_ = amplitude;
		offset_ = offset;
	}

	@Override
	public double computeActivity(double time)
	{
		double val = offset_ + amplitude_ * Math.sin(2*Math.PI*frequency_ * time + phase_);
		return val;
	}

	@Override
	public void compileParameters(ArrayList<String> names, ArrayList<Double> values) 
	{
		names.add("phase");
		names.add("frequency");
		names.add("amplitude");
		names.add("offset");
		
		values.add(phase_);
		values.add(frequency_);
		values.add(amplitude_);
		values.add(offset_);		
	}
	

	@Override
	public void writeSBML(Model model, Species thisSpecies, Species voidSpecies, ArrayList<String> names, ArrayList<Double> values) 
	{
		double startval = getInitialCondition();
		thisSpecies.setInitialAmount(startval);
		
		KineticLaw kineticLaw = new KineticLaw ();
		String math = "\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n";
		math += "\t\t\t<apply><times/>\n";
		math += "\t\t\t\t<cn>6.2831853071795</cn>\n";
		math += "\t\t\t\t<ci>frequency</ci>\n";
		math += "\t\t\t\t<ci>amplitude</ci>\n";
		math += "\t\t\t\t<apply><cos/>\n";
		math += "\t\t\t\t\t<apply><plus/>\n";
		math += "\t\t\t\t\t\t<ci>phase</ci>\n";
		math += "\t\t\t\t\t\t<apply><times/>\n";
		math += "\t\t\t\t\t\t\t<cn>6.2831853071795</cn>\n";
		math += "\t\t\t\t\t\t\t<ci>frequency</ci>\n";
		math += "\t\t\t\t\t\t\t<ci>_time_</ci>\n";
		math += "\t\t\t\t\t\t</apply>\n";
		math += "\t\t\t\t\t</apply>\n";
		math += "\t\t\t\t</apply>\n";
		math += "\t\t\t</apply>\n";
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
		double val = 2*Math.PI * frequency_ * amplitude_ * Math.cos(2*Math.PI*frequency_* time + phase_);
		return val;
	}

	@Override
	public double getInitialCondition()
	{
		return  computeActivity(0);
	}

	@Override
	public void initialization(ArrayList<String> paramNames, ArrayList<Double> paramValues,
			ArrayList<GrnNode> inputGenes) 
	{
		phase_ = paramValues.get(paramNames.indexOf("phase"));
		frequency_ = paramValues.get(paramNames.indexOf("frequency"));
		amplitude_ = paramValues.get(paramNames.indexOf("amplitude"));
		offset_ = paramValues.get(paramNames.indexOf("offset"));
	
		return;
	}	
}