/******************************************************************************
*  Connor Westcott Student					  Developed for CAP5512 (Dr. Wu)
*  Version 1, Feburary 17, 2020
*******************************************************************************/

import java.io.*;
import java.util.*;
import java.text.*;

public class TSP2 extends FitnessFunction{

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/
	public double [][] cityPoints;
	public String dataSet;
	public static Set<Integer> range;
	// A flag set by the parameters for using either Euclidean or Manhatten distance
	public int distanceFunction;

/*******************************************************************************
*                            STATIC VARIABLES                                  *
*******************************************************************************/


/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public TSP2(int distanceFunction) throws java.io.IOException
	{
		name = "Traveling Salesman Problem";
		range = new HashSet<Integer>();
		this.distanceFunction = distanceFunction;
		
		loadCityData();
	}
	
	
	/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/
	// Load the data from a file
	private void loadCityData() throws java.io.IOException
	{
		String line = "";
		BufferedReader cityInput = new BufferedReader(new FileReader(Parameters.dataInputFileName) ); 
		
		this.dataSet = cityInput.readLine().substring(6);
		
		// Init array to store city data
		if (this.dataSet.equals("berlin52"))
			cityPoints = new double [53][2];
		else if (this.dataSet.equals("rl1323"))
			cityPoints = new double [1324][2];
		else
			cityPoints = new double [49][2];
		
		
		for (int i = 0; i < 5; i++)
			cityInput.readLine(); // consume don't care lines
		
		// Put data into array
		do{
			line = cityInput.readLine();
			
			if (line.equals("EOF")) break;
			
			System.out.println(line);
			String indx = line.split(" ")[0];
			String x = line.split(" ")[1];
			String y = line.split(" ")[2];
			
			this.range.add(Integer.parseInt(indx));
			this.cityPoints[Integer.parseInt(indx)][0] = Double.parseDouble(x);
			this.cityPoints[Integer.parseInt(indx)][1] = Double.parseDouble(y);
			
		}while(!line.equals("EOF"));
	}

//  COMPUTE A CHROMOSOME'S RAW FITNESS *************************************

	public void doRawFitness(Chromo X){

		X.rawFitness = 0;
		Set<Integer> chromoRange = new HashSet<Integer>();
		try{
			for (int z=0; z<Parameters.numGenes * Parameters.geneSize; z++){
				int curCity = X.chromo.charAt(z) - '0';
				//System.out.println(X.chromo.charAt(z) - '0');
				//int curCity = Character.getNumericValue(X.chromo.charAt(z));
				//System.out.println(curCity);
				int nextCity = 0;
				if (z+1 < Parameters.numGenes * Parameters.geneSize)
					nextCity = X.chromo.charAt(z+1) - '0';
				else
					nextCity = X.chromo.charAt(0) - '0'; // loop back to the first city
				
				chromoRange.add(curCity); // keep track of the cities we see
				// Distance between 2 points
				if(this.distanceFunction == 0) {
					X.rawFitness += Math.sqrt( 
									Math.pow((this.cityPoints[curCity][0]-this.cityPoints[nextCity][0]), 2)
									+ Math.pow((this.cityPoints[curCity][1]-this.cityPoints[nextCity][1]), 2)
									);
				} else {
					X.rawFitness += Math.abs(this.cityPoints[curCity][0]-this.cityPoints[nextCity][1]) 
									+ Math.abs(this.cityPoints[curCity][0]-this.cityPoints[nextCity][1]);			
				}
			}
			
			Set<Integer> difference = new HashSet<Integer>();
			difference.addAll(TSP2.range);
			difference.removeAll(chromoRange);
			
			if (!difference.isEmpty())
				X.rawFitness = Integer.MAX_VALUE;
		}
		catch(IndexOutOfBoundsException e)
		{
			//System.out.println(e);
			X.rawFitness = Integer.MAX_VALUE;
		}
		
	}

//  PRINT OUT AN INDIVIDUAL GENE TO THE SUMMARY FILE *********************************

	public void doPrintGenes(Chromo X, FileWriter output) throws java.io.IOException{

		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getGeneAlpha(i),11,output);
		}
		output.write("   RawFitness");
		output.write("\n        ");
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getPosIntGeneValue(i),11,output);
		}
		Hwrite.right((int) X.rawFitness,13,output);
		output.write("\n\n");
		return;
	}

/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

}   // End of TSP.java ******************************************************