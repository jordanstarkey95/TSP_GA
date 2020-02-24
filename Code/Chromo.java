/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.util.*;
import java.text.*;

public class Chromo
{
/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	public String chromo;
	public double rawFitness;
	public double sclFitness;
	public double proFitness;
	public int tspRep;
	public int numCities; //use this for random integer range generation

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/

	private static double randnum;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public Chromo(){

		//  Set gene values to a randum sequence of 1's and 0's
		char geneBit;
		chromo = "";
		for (int i=0; i<Parameters.numGenes; i++){
			for (int j=0; j<Parameters.geneSize; j++){
				randnum = Search.r.nextDouble();
				if (randnum > 0.5) geneBit = '0';
				else geneBit = '1';
				this.chromo = chromo + geneBit;
			}
		}

		this.rawFitness = -1;   //  Fitness not yet evaluated
		this.sclFitness = -1;   //  Fitness not yet scaled
		this.proFitness = -1;   //  Fitness not yet proportionalized
	}
	
	public Chromo(boolean isTSP, int numCities, boolean rep2)
	{
		
		chromo = "";
		String gene = "";
		int firstCity = 0;
		Set<Integer> usedCities = new HashSet<Integer>();
		
		
		ArrayList<Integer> randnums = new ArrayList<Integer>();
		
		for (int i = 1; i <= numCities; i++) 
		{
			randnums.add(i);
		}		
		Collections.shuffle(randnums);
		
		for (int location : randnums)
		{
			//this.chromo += Character.toChars(location);
			this.chromo += (char)(location + '0');
			//System.out.println((char)(location + '0'));
			//this.chromo += String.valueOf(location);
			//System.out.println(String.valueOf(location));
		}
		//System.out.println(this.chromo);
		//System.out.println("LENGTH = " + this.chromo.length());
		this.tspRep = 2;
		this.numCities = numCities;
		// set other member variables
		this.rawFitness = -1;   //  Fitness not yet evaluated
		this.sclFitness = -1;   //  Fitness not yet scaled
		this.proFitness = -1;   //  Fitness not yet proportionalized
	}
	
	public Chromo(boolean isTSP, int numCities)
	{			
		chromo = "";
		String gene = "";
		int firstCity = 0;
		Set<Integer> usedCities = new HashSet<Integer>();
		
		
		ArrayList<Integer> randnums = new ArrayList<Integer>();
		
		for (int i = 1; i <= numCities; i++) 
		{
			randnums.add(i);
		}
		
		
		Collections.shuffle(randnums);
		
		
		for (int location : randnums)
		{
			
			if (location <= 15) // 1 Hex Char
			{
				gene = "00" + Integer.toHexString(location);
				this.chromo += gene;
			}
			else if(location <= 255) // 2 Hex Chars
			{
				gene = "0" + Integer.toHexString(location);
				this.chromo += gene;
			}
			else // 3 Hex Chars
			{
				gene = Integer.toHexString(location);
				this.chromo += gene;
			}
		} // end loop
		
		this.tspRep = 1;
		this.numCities = numCities;
		// set other member variables
		this.rawFitness = -1;   //  Fitness not yet evaluated
		this.sclFitness = -1;   //  Fitness not yet scaled
		this.proFitness = -1;   //  Fitness not yet proportionalized
		
	}


/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/

	//  Get Alpha Represenation of a Gene **************************************

	public String getGeneAlpha(int geneID){
		int start = geneID * Parameters.geneSize;
		int end = (geneID+1) * Parameters.geneSize;
		String geneAlpha = this.chromo.substring(start, end);
		return (geneAlpha);
	}

	//  Get Integer Value of a Gene (Positive or Negative, 2's Compliment) ****

	public int getIntGeneValue(int geneID){
		String geneAlpha = "";
		int geneValue;
		char geneSign;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=1; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		geneSign = geneAlpha.charAt(0);
		if (geneSign == '1') geneValue = geneValue - (int)Math.pow(2.0, Parameters.geneSize-1);
		return (geneValue);
	}

	//  Get Integer Value of a Gene (Positive only) ****************************

	public int getPosIntGeneValue(int geneID){
		String geneAlpha = "";
		int geneValue;
		char geneBit;
		geneValue = 0;
		geneAlpha = getGeneAlpha(geneID);
		for (int i=Parameters.geneSize-1; i>=0; i--){
			geneBit = geneAlpha.charAt(i);
			if (geneBit == '1') geneValue = geneValue + (int) Math.pow(2.0, Parameters.geneSize-i-1);
		}
		return (geneValue);
	}

	//  Mutate a Chromosome Based on Mutation Type *****************************

	public void doMutation(){

		String mutChromo = "";
		char x;

		switch (Parameters.mutationType){

		case 1:     //  Replace with new random number

			for (int j=0; j<(Parameters.geneSize * Parameters.numGenes); j++){
				x = this.chromo.charAt(j);
				randnum = Search.r.nextDouble();
				if (randnum < Parameters.mutationRate){
					if (x == '1') x = '0';
					else x = '1';
				}
				mutChromo = mutChromo + x;
			}
			this.chromo = mutChromo;
			break;

		case 2: 	// Replace character for Unicode Representation
			for (int j=0; j<(Parameters.geneSize * Parameters.numGenes); j++){
				x = this.chromo.charAt(j);
				// Generate a random number in the range of the number of cities so that it is valid
				randnum = Search.r.nextInt(numCities);
				if (randnum < Parameters.mutationRate){
					x = (char)(randnum + '0');
				}
				mutChromo = mutChromo + x;
			}
			this.chromo = mutChromo;
			break; 

		default:
			//System.out.println("ERROR - No mutation method selected");
		}
	}

/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

	//  Select a parent for crossover ******************************************

	public static int selectParent(){

		double rWheel = 0;
		int j = 0;
		int k = 0;

		switch (Parameters.selectType){

		case 1:     // Proportional Selection
			randnum = Search.r.nextDouble();
			for (j=0; j<Parameters.popSize; j++){
				rWheel = rWheel + Search.member[j].proFitness;
				if (randnum < rWheel) return(j);
			}
			break;

		case 3:     // Random Selection
			randnum = Search.r.nextDouble();
			j = (int) (randnum * Parameters.popSize);
			return(j);

		case 2:     //  Tournament Selection
			int idx1 = Math.abs(Search.r.nextInt());
			int idx2 = Math.abs(Search.r.nextInt()); 
			double P = 0.5;
			
			if (P > Search.r.nextDouble())
				return (idx1 % Parameters.popSize);
			else
				return (idx2 % Parameters.popSize);
		
		case 4: 	// Rank Selection
			randnum = Search.r.nextDouble();
			
			// Select something
			for (j=0; j<Parameters.popSize; j++){
				rWheel = rWheel + Search.member[j].proFitness;
				//System.out.println(rWheel +","+ randnum);
				if (randnum < rWheel) return(j);
			}
			
			break;

		default:
			//System.out.println("ERROR - No selection method selected");
			break;
		}
	return(-1);
	}

	//  Produce a new child from two parents  **********************************
	private static void validateChildren(Chromo X)
	{
		Set<Integer> chromoRange = new HashSet<Integer>();
		ArrayList<Integer> citiesToFix = new ArrayList<Integer>();
		for (int z=0; z<Parameters.numGenes * Parameters.geneSize; z+=3)
		{
			int curCity = Integer.parseInt(X.chromo.substring(z, z+3), 16);
			
			// Check to make sure the city is in range
			if (curCity > Parameters.numGenes || curCity <= 0)
				citiesToFix.add(z);
			// Check to see if city is repear
			else if (chromoRange.contains(curCity))
				citiesToFix.add(z);
			// Valid Unique city
			else
				chromoRange.add(curCity);			
		}
		Set<Integer> difference = new HashSet<Integer>();
		difference.addAll(TSP.range);
		difference.removeAll(chromoRange);
		if(!difference.isEmpty())
		{
			Iterator<Integer> iter = difference.iterator();
			StringBuffer temp = new StringBuffer(X.chromo);
			String gene = "";
			for (int city : citiesToFix)
			{
				int validCity = iter.next();
				if (validCity <= 15) // 1 Hex Char
				{
					gene = "00" + Integer.toHexString(validCity);
					temp.replace(city, city+3, gene);
				}
				else if(validCity <= 255) // 2 Hex Chars
				{
					gene = "0" + Integer.toHexString(validCity);
					temp.replace(city, city+3, gene);
				}
				else // 3 Hex Chars
				{
					gene = Integer.toHexString(validCity);
					temp.replace(city, city+3, gene);
				}	
			}
			X.chromo = temp.toString();
		}
	}

	private static void validateChildrenTSP2(Chromo X)
	{
		Set<Integer> chromoRange = new HashSet<Integer>();
		ArrayList<Integer> citiesToFix = new ArrayList<Integer>();
		//System.out.println("Before Vvalidate: " + X.chromo);
		//System.out.println("Lenght - " + X.chromo.length());
		for (int z=0; z< Parameters.numGenes; z++)
		{
			int curCity = X.chromo.charAt(z) - '0';
			//System.out.println("Curvalue: " + curCity + " = " + X.chromo.charAt(z));
			
			// Check to make sure the city is in range
			if (curCity > Parameters.numGenes || curCity <= 0)
				citiesToFix.add(z);
			// Check to see if city is repear
			else if (chromoRange.contains(curCity)) {
				citiesToFix.add(z);
			}
			// Valid Unique city
			else
				chromoRange.add(curCity);			
		}
		Set<Integer> difference = new HashSet<Integer>();
		difference.addAll(TSP2.range);
		difference.removeAll(chromoRange);
		//System.out.println("Valid cities: " + chromoRange);
		//System.out.println("Invalid Cities: " + citiesToFix);
		//System.out.println("Remaing valid cities: " + difference);
		if(!difference.isEmpty())
		{
			Iterator<Integer> iter = difference.iterator();
			StringBuffer temp = new StringBuffer(X.chromo);
			String gene = "";
			for (int city : citiesToFix)
			{
				int validCity = iter.next();
				char encodedValidCity = (char)(validCity + '0');
				temp.setCharAt(city, encodedValidCity);
			}
			X.chromo = temp.toString();
		}

		//System.out.println("After validate: " + X.chromo);

		// Test to see that each city is in the validate chromosome
		//for (int z=0; z< Parameters.numGenes; z++)
		//{
		//	int curCity = X.chromo.charAt(z) - '0';
			//System.out.println(curCity);	
		//}
	}

	public static void mateParents(int pnum1, int pnum2, Chromo parent1, Chromo parent2, Chromo child1, Chromo child2){

		int xoverPoint1;
		int xoverPoint2;
		//System.out.println("Parent 1: " + parent1.chromo.length() + " = " + parent1.chromo);
		//System.out.println("Parent 2: " + parent2.chromo.length() + " = " + parent2.chromo);

		switch (Parameters.xoverType){

		case 1:     //  Single Point Crossover

			//  Select crossover point
			xoverPoint1 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));

			//System.out.println(xoverPoint1);

			//  Create child chromosome from parental material
			child1.chromo = parent1.chromo.substring(0,xoverPoint1) + parent2.chromo.substring(xoverPoint1);
			child2.chromo = parent2.chromo.substring(0,xoverPoint1) + parent1.chromo.substring(xoverPoint1);

			//System.out.println("Before Validate: " + child1.chromo);
			//System.out.println("Before Validate: " + child2.chromo);
			
			if (Parameters.problemType.equals("TSP"))
			{
				validateChildren(child1);
				validateChildren(child2);
			}
			else if (Parameters.problemType.equals("TSP2"))
			{
				validateChildrenTSP2(child1);
				validateChildrenTSP2(child2);
			}
			break;

		case 2:     //  Two Point Crossover: Adapted from 1-point crossover implementation.
			
			// Select 2 random crossover points
			xoverPoint1 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));
			xoverPoint2 = 1 + (int)(Search.r.nextDouble() * (Parameters.numGenes * Parameters.geneSize-1));

			// Create children from chromosomes of both parents
			if(xoverPoint1 < xoverPoint2)
			{
			child1.chromo = parent1.chromo.substring(0, xoverPoint1) + parent2.chromo.substring(xoverPoint1, xoverPoint2) + parent1.chromo.substring(xoverPoint2);
			child2.chromo = parent2.chromo.substring(0, xoverPoint1) + parent1.chromo.substring(xoverPoint1, xoverPoint2) + parent2.chromo.substring(xoverPoint2);
			}
			else if(xoverPoint1 >= xoverPoint2)
			{
			child1.chromo = parent1.chromo.substring(0, xoverPoint2) + parent2.chromo.substring(xoverPoint2, xoverPoint1) + parent1.chromo.substring(xoverPoint1);
			child2.chromo = parent2.chromo.substring(0, xoverPoint2) + parent1.chromo.substring(xoverPoint2, xoverPoint1) + parent2.chromo.substring(xoverPoint1);
			}
			else
			{
				System.out.println("ERROR - Bad crossover points selected");
				break;
			}

			// Validate if child chromosomes are valid strings for hex representation
			if (Parameters.problemType.equals("TSP"))
			{
				validateChildren(child1);
				validateChildren(child2);
			}
			//Validate child chromosomes for unicode representation.
			else if (Parameters.problemType.equals("TSP2"))
			{
				validateChildrenTSP2(child1);
				validateChildrenTSP2(child2);
			}
			break;

			
		case 3:     //  Uniform Crossover 

			//System.out.println("Aft before Validate: " + child1.chromo);
			//System.out.println("Afe: " + child2.chromo);

			StringBuffer tempChild1 = new StringBuffer();
			StringBuffer tempChild2 = new StringBuffer();

			for(int i = 0; i < (Parameters.numGenes * Parameters.geneSize); i++)
			{
				randnum = Search.r.nextDouble();
				if(randnum > Parameters.xoverRate)
				{
					//child1.chromo += parent2.chromo.charAt(i);
					//child2.chromo += parent1.chromo.charAt(i);
					tempChild1.append(parent2.chromo.charAt(i));
					tempChild2.append(parent1.chromo.charAt(i));
				}
				else
				{
					//child1.chromo += parent1.chromo.charAt(i);
					//child2.chromo += parent2.chromo.charAt(i);
					tempChild1.append(parent1.chromo.charAt(i));
					tempChild2.append(parent2.chromo.charAt(i));					
				}
			}

			child1.chromo = tempChild1.toString();
			child2.chromo = tempChild2.toString();

			//System.out.println("After uniform and before Validate: " + child1.chromo);
			//System.out.println("After uniform and before Validate: " + child2.chromo);

			// Validate if child chromosomes are valid strings for hex representation
			if (Parameters.problemType.equals("TSP"))
			{
				validateChildren(child1);
				validateChildren(child2);
			}
			//Validate child chromosomes for unicode representation.
			else if (Parameters.problemType.equals("TSP2"))
			{
				validateChildrenTSP2(child1);
				validateChildrenTSP2(child2);
			}
			break;
			
		
		default:
			System.out.println("ERROR - Bad crossover method selected");
		}

		//System.out.println("After Child 1: " + child1.chromo.length() + " + " + child1.chromo);
		//System.out.println("After Child 2: " + child2.chromo.length() + " + " + child2.chromo);

		//  Set fitness values back to zero
		child1.rawFitness = -1;   //  Fitness not yet evaluated
		child1.sclFitness = -1;   //  Fitness not yet scaled
		child1.proFitness = -1;   //  Fitness not yet proportionalized
		child2.rawFitness = -1;   //  Fitness not yet evaluated
		child2.sclFitness = -1;   //  Fitness not yet scaled
		child2.proFitness = -1;   //  Fitness not yet proportionalized
	}

	//  Produce a new child from a single parent  ******************************

	public static void mateParents(int pnum, Chromo parent, Chromo child){

		//  Create child chromosome from parental material
		child.chromo = parent.chromo;

		//  Set fitness values back to zero
		child.rawFitness = -1;   //  Fitness not yet evaluated
		child.sclFitness = -1;   //  Fitness not yet scaled
		child.proFitness = -1;   //  Fitness not yet proportionalized
	}

	//  Copy one chromosome to another  ***************************************

	public static void copyB2A (Chromo targetA, Chromo sourceB){

		targetA.chromo = sourceB.chromo;

		targetA.rawFitness = sourceB.rawFitness;
		targetA.sclFitness = sourceB.sclFitness;
		targetA.proFitness = sourceB.proFitness;
		return;
	}

}   // End of Chromo.java ******************************************************
