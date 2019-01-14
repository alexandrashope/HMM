/**
 * @author: Alexandra Shope
 * @since: June 1, 2018
 * This program uses a Hidden Markov Model and movie scene data to predict desired editing sequences
 * Dr. William Bares contributed to this code. 
 */

import java.util.*;
import javafx.util.Pair;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class BlockingAndShots {
	private static final String COMMA_DELIMITER = ",";
	private static String filePath = ""; //training data file path 
	private static ArrayList<Vector<String>> dataMatrix = new ArrayList<Vector<String>>();
	private static ArrayList<Vector<String>> statesAndObs = new ArrayList<Vector<String>>();
	
	/**  
	 * Initialize matrix to store useful data from files (subject, verb, object, shot type) 
	 **/
	private static void initializeDataMatrix(){
		Vector<String> subjects = new Vector<String>(); //subject, verb, object columns used to construct story event observations
		Vector<String> verb = new Vector<String>();
		Vector<String> object = new Vector<String>();
		Vector<String> standOrSit = new Vector<String>(); //blocking columns used to construct states 
		Vector<String> distance = new  Vector<String>();
		Vector<String> direction = new Vector<String>();
		Vector<String> cameraShots = new Vector<String>();
	    dataMatrix.add(subjects); 
	    dataMatrix.add(verb);
	    dataMatrix.add(object);
	    dataMatrix.add(standOrSit);
	    dataMatrix.add(distance);
	    dataMatrix.add(direction);
	    dataMatrix.add(cameraShots);	    
	}
	
	/**  
	 * Initialize matrix to store states and observations
	 **/
	public static void initializeStatesAndObsMatrix(){
		Vector<String> states = new Vector<String>();
		Vector<String> observations = new Vector<String>();
		statesAndObs.add(states);
		statesAndObs.add(observations);
	}
	
	
	/**  
	 * Parse CSV file and populate data matrix 
	 **/
	public static void parseCSV(String fileName) throws IOException{
		initializeDataMatrix();
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line1 = null; //skip header in file
		br.readLine();

		while ((line1 = br.readLine()) != null) { 
		    String[] cols = line1.split(COMMA_DELIMITER); //parse by column 
		    dataMatrix.get(0).add(cols[0]); //first column: "subjects"
		    dataMatrix.get(1).add(cols[1]); //column 2: "verbs"
		    dataMatrix.get(2).add(cols[2]); //column 3: "objects:
		    dataMatrix.get(3).add(cols[5]); //column 6: "standing or sitting"
		    dataMatrix.get(4).add(cols[6]); //column 7: "Distance" 
		    dataMatrix.get(5).add(cols[7]); //column 8: "Direction" 
		    dataMatrix.get(6).add(cols[8]); //column 9: "Camera shots" 		    
		}	
		br.close();
	}
	
	/**  
	 * Loop through test file of file names call to ParseCSV 
	 **/
	public static void parseMultipleCSV(String fileNames){
		ArrayList<String> files = getFileName(fileNames);
		for(int i = 0; i < files.size(); i++){
			try {
				parseCSV(files.get(i));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	
	/**  
	 * Read file containing names of files to parse
	 * @return ArrayList of file names
	 **/
	public static ArrayList<String> getFileName(String fileNames) {
		BufferedReader br;
		ArrayList<String> names = new ArrayList<String>(); //store names of files 
		String line;
		try {
			br = new BufferedReader(new FileReader(fileNames));
			while ((line = br.readLine()) != null) {
	            names.add(line); 
	        }
			
			br.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return names;
	}
		
	
	/**  
	 * Parse test file to construct generic observations 
	 **/
	public static void parseCSVTestObservations(String fileName){
		BufferedReader br;
		String line1 = null; //skip header in file
		dataMatrix.clear();
		initializeDataMatrix(); 
		try {
			br = new BufferedReader(new FileReader(fileName));
			br.readLine();
			while ((line1 = br.readLine()) != null) { 
			    String[] cols = line1.split(COMMA_DELIMITER); //parse by column 
			    dataMatrix.get(0).add(cols[0]);
			    dataMatrix.get(1).add(cols[1]);
			    dataMatrix.get(2).add(cols[2]);
			}
			br.close();
			
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**  
	 * @return sequence of movie scene subjects  
	 **/
	public static  Vector<String> getSubjectSequence(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> subjects = data.get(0); 
		return subjects;		
	}	
	
	/**  
	 * @return sequence of movie scene verbs  
	 **/
	public static Vector<String> getVerbsSequence(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> verbs = data.get(1);
		return verbs;		
	}
	
	/**  
	 * @return sequence of movie scene objects  
	 **/
	public static Vector<String> getObjectsSequence(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> objects = data.get(2);
		return objects;		
	}
	
	/**  
	 * @return sequence of movie scene objects  
	 **/
	public static Vector<String> getStandOrSit(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> standOrSit = data.get(3); 
		return standOrSit;
	}
	
	/**  
	 * @return sequence of blocking distances 
	 **/
	public static Vector<String> getDistance(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> distance = data.get(4); 
		return distance;
	}
	
	/**  
	 * @return sequence of blocking directions
	 **/
	public static Vector<String> getDirection(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> direction = data.get(5); 
		return direction;
	}
	
	/**  
	 * @return sequence of camera shots
	 **/
	public static Vector<String> getCameraShots(){
		ArrayList<Vector<String>> data = dataMatrix;
		Vector<String> camShots = data.get(6); 
		return camShots;
	}
	
	
	/**  
	 * @set training observations as story events to predict blocking
	 **/
	public static void setObservationsSequenceBlocking(){
		statesAndObs.get(1).clear(); //observations vector in statesAndObs Matrix cleared
		Vector<String> subjects = getSubjectSequence();
		Vector<String> verbs = getVerbsSequence();
		Vector<String> objects = getObjectsSequence();
		
		Vector<String> observations = new Vector<String>();
		for(int i = 0; i < getSubjectSequence().size(); i++){ //subject sequence is the same size as verbs and objects 
			observations.add(subjects.get(i) + " " + verbs.get(i) + " " + objects.get(i)); //construct observations
		}
		statesAndObs.get(1).addAll(observations); //set observations
	}
	
	
	/**  
	 * @set training observations as story events and corresponding blocking to predict camera shots
	 **/
	public static void setObservationsCameraShots(){
		statesAndObs.get(1).clear();  //observations vector in statesAndObs Matrix cleared
		Vector<String> subjects = getSubjectSequence();
		Vector<String> verbs = getVerbsSequence();
		Vector<String> objects = getObjectsSequence();
		Vector<String> standSit = getStandOrSit();
		Vector<String> distance = getDistance();
		Vector<String> direction = getDirection();
		
		Vector<String> obs = new Vector<String>();
		for(int i = 0; i < getSubjectSequence().size(); i++){ 
			obs.add(subjects.get(i) + " " + verbs.get(i) + " " + objects.get(i)
			+ " " + standSit.get(i) + "-" + distance.get(i) + "-" + direction.get(i)); //construct observations
		}
		statesAndObs.get(1).addAll(obs);//set observations
	}
	
	/**  
	 * @return observations sequence
	 **/
	public static Vector<String> getObservationSequence(){
		return statesAndObs.get(1);
	}
	
	/**  
	 * @set states sequence as blocking for first HMM 
	 **/
	public static void setStatesSequenceBlocking(){
		statesAndObs.get(0).clear(); //states vector in statesAndObs Matrix cleared
		
		Vector<String> standOrSit = getStandOrSit(); //data to construct blocking states 
		Vector<String> distance = getDistance();
		Vector<String> direction = getDirection();
		Vector<String> stateSeq = new Vector<String>();
		
		for(int i = 0; i < standOrSit.size(); i++){
			stateSeq.add(standOrSit.get(i)+ "-" + distance.get(i) + "-" + direction.get(i)); //construct states
		} 	
		statesAndObs.get(0).addAll(stateSeq); //set states
	}
	
	/**  
	 * @set states sequence as camera shots for second HMM 
	 **/
	public static void setStatesSequenceCameraShots(){
		statesAndObs.get(0).clear(); //states vector in statesAndObs Matrix cleared
		
		Vector<String> camShots = getCameraShots();
		statesAndObs.get(0).addAll(camShots); //set states
	}
	
	/**  
	 * @return states sequence 
	 **/
	public static Vector<String> getStatesSequence(){
		Vector<String> stateSeq = statesAndObs.get(0);
		return stateSeq;
	}
	
	
	/**  
	 * @return possible hidden states 
	 **/
	public static Vector<String> getStates(){
		Vector<String> statesSequence = getStatesSequence(); //current states sequence 
		
		Vector<String> states = new Vector<String>(); //vector of unique states
		for(int i=0; i < statesSequence.size(); i++){
			if(!states.contains(statesSequence.get(i))){ //if state has not been added to new vector, add it, else skip
				states.add(statesSequence.get(i));		
			}
		}
		return states;
	}
	
	
	/**  
	 * Calculates initial probabilities; assuming all equally likely
	 * @return probabilities of each shot being the initial shot
	 **/
	public static Hashtable<String,Double> calculateInitialProbabilities(){
		Hashtable<String,Double> initialProbabilities = new Hashtable<String,Double>();
		Vector<String> states = getStates();
		double numStates = (double) states.size();
		double prob = 1.0 / numStates;
		for(int i = 0; i < numStates; i++){
			initialProbabilities.put(states.get(i), prob);
		}
		return initialProbabilities;
	}		

	/**  
	 * Calculate emission probabilities and populate hashtable
	 * @return emissionMatrix hashtable
	 **/ 
	public static Hashtable<Pair<String,String>,Double> calculateEmissionProbs(){
		Hashtable<Pair<String, String>, Double> emissionMatrix = new Hashtable<Pair<String, String>, Double>();	 
		Vector<String> states = getStates();
		Vector<String> statesSeq = getStatesSequence();		
		Vector<String> obsSeq = getObservationSequence();
				
		Vector<Integer> indeces = new Vector<Integer>(); //vector that stores indeces each time hidden state occurs in states sequence
		double prob = 0.0; //emission probability 
		
		ArrayList<String> correspondingSeq; //sequence of corresponding observations for each state 
		Pair<String,String> emPair = new Pair<String,String>(null, null);
		Pair<String,String> emPair2 = new Pair<String,String>(null, null);		
		
		int statesCounter = 0; //count number of occurrences of each state in state sequence 
		int obsCounter; //count number of times corresponding observations occur
		
		for(int i = 0; i < states.size(); i++){ //loop through possible states
			indeces.clear(); 
			for(int j = 0; j < statesSeq.size(); j++){ //loop through states sequence 
				if(states.get(i).compareTo(statesSeq.get(j)) == 0){ //find indeces in sequence where state occurs and store them
					indeces.add(j); 
				}
	 	
			} //END j
			statesCounter = indeces.size();	
			
			correspondingSeq = new ArrayList<String>(indeces.size()); //sequence of corresponding observations 
			for(int k = 0; k < indeces.size(); k++){
				correspondingSeq.add(obsSeq.get(indeces.get(k))); //store corresponding observations
			} //END k
			
			for(int l = 0; l < correspondingSeq.size(); l++){
				obsCounter = Collections.frequency(correspondingSeq,correspondingSeq.get(l)); //count # of times observation occurs given hidden state 
				prob = (double) obsCounter / (double) statesCounter; //calculate probability				
				emPair = new Pair<String,String>(states.get(i),correspondingSeq.get(l)); //construct emission pair 
				emissionMatrix.put(emPair, prob); //populate hashtable 
			} //END l
			
			for(int m = 0; m < obsSeq.size(); m++){
				if(!correspondingSeq.contains(obsSeq.get(m))){ 
					emPair2 = new Pair<String,String>(states.get(i), obsSeq.get(m)); 
					emissionMatrix.put(emPair2, 0.0); //populate hashtable with pairs that did not occur as probability 0%
				}
			} // END m		
		} // END i
		return emissionMatrix;
	} 
	
	/**  
	 * Helper method to populate transition matrix
	 * @param hidden state
	 * @return index of hidden state in vector of possible states
	 **/ 
	public static int findStateIndex(String state){
		Vector<String> states = getStates();
		return states.indexOf(state);
	}
	
	/**  
	 * Calculate transition probabilities and populate hashtable
	 * @return transitionMatrix hashtable
	 **/ 
	public static Hashtable<Pair<String,String>,Double> calculateTransitionProbs(){
		Hashtable<Pair<String, String>, Double> transitionMatrix = new Hashtable<Pair<String, String>, Double>();
		Vector<String> states = getStates();
		Vector<String> statesSeq = getStatesSequence();
			
		double[][] transitions = new double[states.size()][states.size()]; //transition matrix to store number of transition between pairs
		String fromState;
		String toState;
		int fromStateIndex; //to update matrix;
		int toStateIndex; //to update matrix; 
		int numTransitions = 0;//this number is wrong, fix this!!!
		
		//count transitions, update matrix
		for(int i = 0; i < statesSeq.size()-1; i++){
			fromState = statesSeq.get(i);
			toState = statesSeq.get(i+1);
			fromStateIndex = findStateIndex(fromState);
			toStateIndex = findStateIndex(toState);
			transitions[fromStateIndex][toStateIndex] ++;
		}
		
		//calculate probabilities, update matrix 
		for(int i = 0; i < states.size(); i++){
			numTransitions = 0;
			for(int j = 0; j < states.size(); j++){
				
				if(transitions[i][j] > 0.0){
					numTransitions++;
				}	

			} //END J loop 
			for(int j = 0; j < states.size(); j++){
				transitions[i][j] /= (double) numTransitions;
			} // J loop 
		} //END loop 
		
		//populate hashtable 
		for(int i = 0; i < transitions[0].length; i++){ 
			for(int j = 0; j < transitions.length; j++){
				transitionMatrix.put(new Pair<String,String>(states.get(i),states.get(j)),transitions[i][j]);
			}
		} //END loop 	
		return transitionMatrix;
	}
	
	/**  
	 * Build and train Hidden Markov Model to predict blocking
	 * @return HiddenMarkovModel 
	 **/ 
	public static HiddenMarkovModel trainBlocking(){	
		parseMultipleCSV(filePath); //data matrix populated, use this to populate states and observations matrix
		initializeStatesAndObsMatrix(); //initialize states and observations matrix 
		setStatesSequenceBlocking(); //pick states as blocking out of data matrix and store states sequence 
		setObservationsSequenceBlocking(); //pick story events as observations out of data matrix
		
		Vector<String> states = getStates(); //unique states 
		Vector<String> observationSequence = getObservationSequence(); //observation sequence 
		
		Hashtable<String,Double> initialProbabilities = calculateInitialProbabilities();
		Hashtable<Pair<String,String>,Double> emissionMatrix = calculateEmissionProbs();
		Hashtable<Pair<String,String>,Double> transitionMatrix = calculateTransitionProbs();
		
		HiddenMarkovModel trained_Blocking = new HiddenMarkovModel("Blocking",states, observationSequence,initialProbabilities,transitionMatrix,emissionMatrix);
		return trained_Blocking;
	}
	
	
	/**  
	 * Predict and display suggested blocking sequence 
	 * @return newObservations 
	 **/ 
	public static Vector<String> blockingHMMPredict(String predictPath){
		HiddenMarkovModel hmm_Blocking = trainBlocking();
		Scanner keyboard = new Scanner(System.in);
		
		Vector<String> states = getStates(); //blocking states
		dataMatrix.clear();
		initializeDataMatrix();
		statesAndObs.clear();
		initializeStatesAndObsMatrix();
		
		parseCSVTestObservations(predictPath);
		setObservationsSequenceBlocking();
		Vector<String> observations = getObservationSequence();
		
		Vector<String> blocking = new Vector<String>(); //used to store predicted blockings 
		String path = hmm_Blocking.getOptimalStateSequenceUsingViterbiAlgorithm(states, observations);
		System.out.println("\nYour story events: " + observations);
		System.out.println("Suggested blocking sequence: " + path);
		System.out.println("\nAre you statisfied with this blocking sequence?\n"
				+ " If so, type Y and a sequence of camera shots will automatically be computed.\n If not, type N to edit the sequence.");
		
		String answer = keyboard.nextLine();
		if(answer.compareTo("Y") == 0 || answer.compareTo("y") == 0){
			String[] blockingList = path.split(" -> ");
			for(int i = 0; i < blockingList.length -1; i++){ //ignore total cost at end 
				blocking.add(blockingList[i]);
			}
		}

		else if(answer.compareTo("N") == 0 || answer.compareTo("n") == 0){
			System.out.println("Type your edited sequence in the same format as above (case sensitive),"
					+ " separated by commas. For example: standing-close-toward, standing-medium-toward, standing-far-away");
			String newEdits = keyboard.nextLine();
			String[] blockingList = newEdits.split(", ");
			for(int i = 0; i < blockingList.length; i++){ //ignore total cost at end 
				blocking.add(blockingList[i]);
			}
		}
		Vector<String> newObservations = new Vector<String>();
		for(int i = 0; i < observations.size(); i++){
			newObservations.add(observations.get(i) + " " + blocking.get(i));
		}
		keyboard.close();
		return newObservations;
	}	
	
	/**  
	 * Helper method to retrieve new observations for prediction
	 * @return newTestObs
	 **/ 
	public static Vector<String> getFirstHMMData(String predictPath){
		Vector<String> newTestObs = blockingHMMPredict(predictPath);
		return newTestObs; 
	}
	
	/**  
	 * Build and train Hidden Markov Model to predict Camera Shots
	 * @return HiddenMarkovModel 
	 **/ 
	public static HiddenMarkovModel trainCameraShots(){
		dataMatrix.clear();
		initializeDataMatrix();
		statesAndObs.clear();
		initializeStatesAndObsMatrix();
		parseMultipleCSV(filePath); //training on same data
		
		setStatesSequenceCameraShots(); //states are now camera shots 
		setObservationsCameraShots(); //observations are now story events and blocking 
		
		Vector<String> observationSequence = getObservationSequence();
		//System.out.println("training obs: " + observationSequence);
		Vector<String> states = getStates();
		//System.out.println("training states: " + states);

		Hashtable<String,Double> initialProbabilities = calculateInitialProbabilities();
		Hashtable<Pair<String,String>,Double> emissionMatrix = calculateEmissionProbs();
		Hashtable<Pair<String,String>,Double> transitionMatrix = calculateTransitionProbs();
		
		HiddenMarkovModel trained = new HiddenMarkovModel("Camera Shots", states, observationSequence,initialProbabilities,transitionMatrix,emissionMatrix);
		return trained;
	}
	

	/**  
	 * Predict and display sequence of camera shots
	 **/ 
	public static void cameraShotsHMMPredict(String predictPath){
		Vector<String> observationSequence = getFirstHMMData(predictPath); //observations to predict on 
		HiddenMarkovModel hmm_CameraShots = trainCameraShots();
		Vector<String> states = hmm_CameraShots.getStates();
		String path = hmm_CameraShots.getOptimalStateSequenceUsingViterbiAlgorithm(states, observationSequence);
		System.out.println("Suggested camera shot sequence: " + path);		
	}
	
		
	public static void main(String [] args){	
		System.out.println("This program reads in a text file containing file paths for training Data.");		
		System.out.println("Enter path of file containing blocking traning data file names: ");
		Scanner keyboard = new Scanner(System.in);
		filePath = keyboard.nextLine();
		System.out.println("Enter the observations file path: ");
		String predictPath = keyboard.nextLine();
		cameraShotsHMMPredict(predictPath);
		keyboard.close();     
	}
}
	

