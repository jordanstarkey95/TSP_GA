/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.util.*;
import java.text.*;

public class Search {

/*******************************************************************************
*                           INSTANCE VARIABLES                                 *
*******************************************************************************/

/*******************************************************************************
*                           STATIC VARIABLES                                   *
*******************************************************************************/

	public static FitnessFunction problem;

	public static Chromo[] member;
	public static Chromo[] child;
	public static Chromo[] sorted;

	public static Chromo bestOfGenChromo;
	public static int bestOfGenR;
	public static int bestOfGenG;
	public static Chromo bestOfRunChromo;
	public static int bestOfRunR;
	public static int bestOfRunG;
	public static Chromo bestOverAllChromo;
	public static int bestOverAllR;
	public static int bestOverAllG;
 
	public static int optimalChromo; // Stores the generation that the best indv was found
	public static int optimalSolution;
	
	public static boolean isTSP;
	public static int tspRep;

	public static double sumRawFitness;
	public static double sumRawFitness2;	// sum of squares of fitness
	public static double sumSclFitness;
	public static double sumProFitness;
	public static double defaultBest;
	public static double defaultWorst;

	public static double averageRawFitness;
	public static double stdevRawFitness;
	public static double stdevAvgAvgFit;

	public static int G;
	public static int R;
	public static Random r = new Random();
	private static double randnum;

	private static int memberIndex[];
	private static double memberFitness[];
	private static int TmemberIndex;
	private static double TmemberFitness;

	private static double fitnessStats[][];  // 0=Avg, 1=Best

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/


/*******************************************************************************
*                             MEMBER METHODS                                   *
*******************************************************************************/


/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

	public static void main(String[] args) throws java.io.IOException{

		Calendar dateAndTime = Calendar.getInstance(); 
		Date startTime = dateAndTime.getTime();

	//  Read Parameter File
		System.out.println("\nParameter File Name is: " + args[0] + "\n");
		Parameters parmValues = new Parameters(args[0]);

	//  Write Parameters To Summary Output File
		String summaryFileName = Parameters.expID + "_summary.txt";
		FileWriter summaryOutput = new FileWriter(summaryFileName);
		parmValues.outputParameters(summaryOutput);

	//	Set up Fitness Statistics matrix
		fitnessStats = new double[2][Parameters.generations];
		for (int i=0; i<Parameters.generations; i++){
			fitnessStats[0][i] = 0;
			fitnessStats[1][i] = 0;
		}

	//	Problem Specific Setup - For new new fitness function problems, create
	//	the appropriate class file (extending FitnessFunction.java) and add
	//	an else_if block below to instantiate the problem.
 
		if (Parameters.problemType.equals("NM")){
				problem = new NumberMatch();
				optimalSolution = 0; // if the numbers are "guessed" correctly then the ideal fitness would be zero based on the function
		}
		else if (Parameters.problemType.equals("OM")){
				problem = new OneMax();
				optimalSolution = Parameters.geneSize; // opt solution is all 1 so fitness would == numGenes for op solution
		}
		else if (Parameters.problemType.equals("TSP"))
		{
			problem = new TSP();
			tspRep = 1;
			isTSP = true;
			optimalSolution = 0;
		}
		else if (Parameters.problemType.equals("TSP2"))
		{
			problem = new TSP2();
			tspRep = 2;
			isTSP = true;
			optimalSolution = 0;
		}
		else System.out.println("Invalid Problem Type");

		System.out.println(problem.name);

	//	Initialize RNG, array sizes and other objects
		r.setSeed(Parameters.seed);
		memberIndex = new int[Parameters.popSize];
		memberFitness = new double[Parameters.popSize];
		member = new Chromo[Parameters.popSize];
		child = new Chromo[Parameters.popSize];
		
		if (isTSP){
			if (tspRep == 1)
			{
				bestOfGenChromo = new Chromo(isTSP, Parameters.numGenes);
				bestOfRunChromo = new Chromo(isTSP, Parameters.numGenes);
				bestOverAllChromo = new Chromo(isTSP, Parameters.numGenes);
			}
			else
			{
				bestOfGenChromo = new Chromo(isTSP, Parameters.numGenes, true);
				bestOfRunChromo = new Chromo(isTSP, Parameters.numGenes, true);
				bestOverAllChromo = new Chromo(isTSP, Parameters.numGenes, true);
			}
		}
		else{
			bestOfGenChromo = new Chromo();
			bestOfRunChromo = new Chromo();
			bestOverAllChromo = new Chromo();
		}
	

		if (Parameters.minORmax.equals("max")){
			defaultBest = 0;
			defaultWorst = 999999999999999999999.0;
			optimalChromo = -1;
		}
		else{
			defaultBest = 999999999999999999999.0;
			defaultWorst = 0;
			optimalChromo = -1;
		}

		bestOverAllChromo.rawFitness = defaultBest;

		//  Start program for multiple runs
		for (R = 1; R <= Parameters.numRuns; R++){
			// Reset vaules to default
			bestOfRunChromo.rawFitness = defaultBest;
			optimalChromo = -1;
			System.out.println("Run: " + R);
			System.out.println();

			//	Initialize First Generation
			for (int i=0; i<Parameters.popSize; i++){
				if (isTSP){
					if (tspRep == 1)
					{
						member[i] = new Chromo(isTSP, Parameters.numGenes);
						child[i] = new Chromo(isTSP, Parameters.numGenes);
					}
					else
					{
						member[i] = new Chromo(isTSP, Parameters.numGenes, true);
						child[i] = new Chromo(isTSP, Parameters.numGenes, true);
					}
				}
				else{
					member[i] = new Chromo();
					child[i] = new Chromo();
				}
			}

			//	Begin Each Run
			for (G=0; G<Parameters.generations; G++){

				sumProFitness = 0;
				sumSclFitness = 0;
				sumRawFitness = 0;
				sumRawFitness2 = 0;
				bestOfGenChromo.rawFitness = defaultBest;

				//	Test Fitness of Each Member
				for (int i=0; i<Parameters.popSize; i++){

					member[i].rawFitness = 0;
					member[i].sclFitness = 0;
					member[i].proFitness = 0;

					problem.doRawFitness(member[i]);

					sumRawFitness = sumRawFitness + member[i].rawFitness;
					sumRawFitness2 = sumRawFitness2 +
						member[i].rawFitness * member[i].rawFitness;

					if (Parameters.minORmax.equals("max")){
						if (member[i].rawFitness > bestOfGenChromo.rawFitness){
							Chromo.copyB2A(bestOfGenChromo, member[i]);
							bestOfGenR = R;
							bestOfGenG = G;
						}
						if (member[i].rawFitness > bestOfRunChromo.rawFitness){
							Chromo.copyB2A(bestOfRunChromo, member[i]);
							bestOfRunR = R;
							bestOfRunG = G;
              if (optimalChromo < 0 && member[i].rawFitness == optimalSolution)
                optimalChromo = G;
						}
						if (member[i].rawFitness > bestOverAllChromo.rawFitness){
							Chromo.copyB2A(bestOverAllChromo, member[i]);
							bestOverAllR = R;
							bestOverAllG = G;
						}
					}
					else {
						if (member[i].rawFitness < bestOfGenChromo.rawFitness){
							Chromo.copyB2A(bestOfGenChromo, member[i]);
							bestOfGenR = R;
							bestOfGenG = G;
						}
						if (member[i].rawFitness < bestOfRunChromo.rawFitness){
							Chromo.copyB2A(bestOfRunChromo, member[i]);
							bestOfRunR = R;
							bestOfRunG = G;
						if (optimalChromo == -1 && member[i].rawFitness == optimalSolution)
							optimalChromo = G;
						}
						if (member[i].rawFitness < bestOverAllChromo.rawFitness){
							Chromo.copyB2A(bestOverAllChromo, member[i]);
							bestOverAllR = R;
							bestOverAllG = G;
						}
					}
				}

				// Accumulate fitness statistics
				fitnessStats[0][G] += sumRawFitness / Parameters.popSize;
				fitnessStats[1][G] += bestOfGenChromo.rawFitness;

				averageRawFitness = sumRawFitness / Parameters.popSize;
				stdevRawFitness = Math.sqrt(
							Math.abs(sumRawFitness2 - 
							sumRawFitness*sumRawFitness/Parameters.popSize)
							/
							(Parameters.popSize-1)
							);

				// Output generation statistics to screen
				System.out.println(R + "\t" + G +  "\t" + (int)bestOfGenChromo.rawFitness + "\t" + averageRawFitness + "\t" + stdevRawFitness);

				// Output generation statistics to summary file
				summaryOutput.write(" R ");
				Hwrite.right(R, 3, summaryOutput);
				summaryOutput.write(" G ");
				Hwrite.right(G, 3, summaryOutput);
				Hwrite.right((int)bestOfGenChromo.rawFitness, 7, summaryOutput);
				Hwrite.right(averageRawFitness, 11, 3, summaryOutput);
				Hwrite.right(stdevRawFitness, 11, 3, summaryOutput);
				summaryOutput.write("\n");


		// *********************************************************************
		// **************** SCALE FITNESS OF EACH MEMBER AND SUM ***************
		// *********************************************************************
		// *********************************************************************
		// ****************    INDV SELECTION ALGORITHMS PREP   ***************
		// *********************************************************************

				switch(Parameters.scaleType){

				case 0:     // No change to raw fitness
					for (int i=0; i<Parameters.popSize; i++){
						member[i].sclFitness = member[i].rawFitness + .000001;
						sumSclFitness += member[i].sclFitness;
					}
					break;

				case 1:     // Fitness not scaled.  Only inverted.
					for (int i=0; i<Parameters.popSize; i++){
						member[i].sclFitness = 1/(member[i].rawFitness + .000001);
						sumSclFitness += member[i].sclFitness;
					}
					break;

				case 2:     // Fitness scaled by Rank (Maximizing fitness)

					//  Copy genetic data to temp array
					for (int i=0; i<Parameters.popSize; i++){
						memberIndex[i] = i;
						memberFitness[i] = member[i].rawFitness;
					}
					//  Bubble Sort the array by floating point number
					for (int i=Parameters.popSize-1; i>0; i--){
						for (int j=0; j<i; j++){
							if (memberFitness[j] > memberFitness[j+1]){
								TmemberIndex = memberIndex[j];
								TmemberFitness = memberFitness[j];
								memberIndex[j] = memberIndex[j+1];
								memberFitness[j] = memberFitness[j+1];
								memberIndex[j+1] = TmemberIndex;
								memberFitness[j+1] = TmemberFitness;
							}
						}
					}
					//  Copy ordered array to scale fitness fields
					for (int i=0; i<Parameters.popSize; i++){
						member[memberIndex[i]].sclFitness = i;
						sumSclFitness += member[memberIndex[i]].sclFitness;
					}

					break;

				case 3:     // Fitness scaled by Rank (minimizing fitness)

					//  Copy genetic data to temp array
					for (int i=0; i<Parameters.popSize; i++){
						memberIndex[i] = i;
						memberFitness[i] = member[i].rawFitness;
					}
					//  Bubble Sort the array by floating point number
					for (int i=1; i<Parameters.popSize; i++){
						for (int j=(Parameters.popSize - 1); j>=i; j--){
							if (memberFitness[j-i] < memberFitness[j]){
								TmemberIndex = memberIndex[j-1];
								TmemberFitness = memberFitness[j-1];
								memberIndex[j-1] = memberIndex[j];
								memberFitness[j-1] = memberFitness[j];
								memberIndex[j] = TmemberIndex;
								memberFitness[j] = TmemberFitness;
							}
						}
					}
					//  Copy array order to scale fitness fields
					for (int i=0; i<Parameters.popSize; i++){
						member[memberIndex[i]].sclFitness = i;
						sumSclFitness += member[memberIndex[i]].sclFitness;
					}

					break;

				default:
					System.out.println("ERROR - No scaling method selected");
				}


		// *********************************************************************
		// ****** PROPORTIONALIZE SCALED FITNESS FOR EACH MEMBER AND SUM *******
		// *********************************************************************

				for (int i=0; i<Parameters.popSize; i++){
					member[i].proFitness = member[i].sclFitness/sumSclFitness;
					sumProFitness = sumProFitness + member[i].proFitness;
				}

		// *********************************************************************
		// ************ CROSSOVER AND CREATE NEXT GENERATION *******************
		// *********************************************************************

				int parent1 = -1;
				int parent2 = -1;

				if (Parameters.selectType == 4)
				{
					System.out.println("Sorting....");
					// Sort in non-decreasing order
					Arrays.sort(member, new Comparator<Chromo>() {
						public int compare(Chromo c1, Chromo c2){
							return Double.valueOf(c1.rawFitness).compareTo(c2.rawFitness);
						}
					});
					System.out.println("Sorting Complete!");
					// get the total rank sum
					int r_sum = 0;
					for (int j = 0; j < Parameters.popSize; j++)
						r_sum += j;
			
					System.out.println("R Sum found...");
			
					// Find the P_s Rank
					for (int j=0; j<Parameters.popSize; j++){
						Search.member[j].proFitness = ((double)j+1)/r_sum; // Rank which is equal to its position in the array (plus 1)
					}
					System.out.println("P_s Found...");
				}
				
				
				//  Assumes always two offspring per mating
				for (int i=0; i<Parameters.popSize; i=i+2){

					//	Select Two Parents
					parent1 = Chromo.selectParent();
					parent2 = parent1;
					while (parent2 == parent1){
						parent2 = Chromo.selectParent();
						//System.out.println(parent1 + "," + parent2);
					}

					//	Crossover Two Parents to Create Two Children
					randnum = r.nextDouble();
					if (randnum < Parameters.xoverRate){
						Chromo.mateParents(parent1, parent2, member[parent1], member[parent2], child[i], child[i+1]);
					}
					else {
						Chromo.mateParents(parent1, member[parent1], child[i]);
						Chromo.mateParents(parent2, member[parent2], child[i+1]);
					}
				} // End Crossover

				//	Mutate Children
				for (int i=0; i<Parameters.popSize; i++){
					child[i].doMutation();
				}

				//	Swap Children with Last Generation
				for (int i=0; i<Parameters.popSize; i++){
					Chromo.copyB2A(member[i], child[i]);
				}

			} //  Repeat the above loop for each generation

			Hwrite.left(bestOfRunR, 4, summaryOutput);
			Hwrite.right(bestOfRunG, 4, summaryOutput);

			problem.doPrintGenes(bestOfRunChromo, summaryOutput);
			Hwrite.left("O", 8, summaryOutput);
			Hwrite.center(optimalChromo, 8, summaryOutput);
			Hwrite.right("\n\n", 3, summaryOutput);

			System.out.println(R + "\t" + "B" + "\t"+ (int)bestOfRunChromo.rawFitness);

		} //End of a Run

		Hwrite.left("B", 8, summaryOutput);

		problem.doPrintGenes(bestOverAllChromo, summaryOutput);
   

		//	Output Fitness Statistics matrix
		summaryOutput.write("Gen                 AvgFit              BestFit \n");
		for (int i=0; i<Parameters.generations; i++){
			Hwrite.left(i, 15, summaryOutput);
			double avgAvgFit = fitnessStats[0][i]/Parameters.numRuns;
			double avgBestFit = fitnessStats[1][i]/Parameters.numRuns;
			
			
			stdevAvgAvgFit= Math.sqrt(
							Math.abs(sumRawFitness2 - 
							sumRawFitness*sumRawFitness/Parameters.popSize)
							/
							(Parameters.popSize-1)
							);
			
			Hwrite.left(avgAvgFit, 20, 2, summaryOutput);
			Hwrite.left(avgBestFit, 20, 2, summaryOutput);
			summaryOutput.write("\n");
		}

		summaryOutput.write("\n");
		summaryOutput.close();

		System.out.println();
		System.out.println("Start:  " + startTime);
		dateAndTime = Calendar.getInstance(); 
		Date endTime = dateAndTime.getTime();
		System.out.println("End  :  " + endTime);

	} // End of Main Class

}   // End of Search.Java ******************************************************

