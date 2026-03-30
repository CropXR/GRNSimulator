package ch.epfl.lis.gnw;

public class GrnNodeFactory 
{
	static GrnNode createNode(String id, GeneNetwork grn)
	{
		GrnNode node = null;
		if (id.contains("ConstantEnv"))
		{
			node = new ConstantEnvironment(grn);
		}
		else if (id.contains("LinearEnv"))
		{
			node = new LinearEnvironment(grn);
		}
		else if (id.contains("OscillatingEnv"))
		{
			node = new OscillatingEnvironment(grn);
		}
		else if (!id.equals("_void_") &&!id.equals("_time_") )
		{
			node = new HillGene(grn);
		}

		return node;
	}
}
