import java.util.*;
import javafx.util.Pair;

public class HiddenMarkovModel {
    private String name;
    private int numberOfStates;
    private int numberOfObservations;
    private Vector<String> states; 
    private Vector<String> observations;
    private Hashtable<String, Double> initialProbabilities;
    private Hashtable<Pair<String, String>, Double> transitionMatrix;
    private Hashtable<Pair<String, String>, Double> emissionMatrix;
	private Pair myPair;

    /**
     * A constructor that initialize the class attributes
     * @param states A Vector that is the states of the model
     * @param observations  A Vector that is the observations of the model
     * @param initialProbabilities A Hashtable that is the initial probability vector of the states
     * @param transitionMatrix A Hashtable the transition matrix between the states
     * @param emissionMatrix A Hashtable that is the emission matrix between the states and the observations
     */

    
    
    public HiddenMarkovModel(String name, Vector<String> states, Vector<String> observations, Hashtable<String, Double> initialProbabilities, Hashtable<Pair<String, String>, Double> transitionMatrix, Hashtable<Pair<String, String>, Double> emissionMatrix){
        this.name = name;
        this.states = states;
        this.numberOfStates = states.size();
        this.observations = observations;
        this.numberOfObservations = observations.size();
        this.initialProbabilities = initialProbabilities;
        this.transitionMatrix = transitionMatrix;
        this.emissionMatrix = emissionMatrix;
    }

    public String getName() {
        return this.name;
    }

    /**
     * Get the number of states in the model
     * @return An integer that specifies the number of states in the model
     */
    public int getNumberOfStates() {
        return this.numberOfStates;
    }

    /**
     * Get the model states
     * @return A Vector which is the states of the model
     */
    public Vector<String> getStates() {
        return states;
    }

    /**
     * Set the number of states in the model
     * @param numberOfStates integer
     */
    public void setNumberOfStates(int numberOfStates) {
        this.numberOfStates = numberOfStates;
    }

    /**
     * Get the number of observations in the model
     * @return An integer that specifies the number of observations in the model
     */
    public int getNumberOfObservations() {
        return numberOfObservations;
    }

    /**
     * Get the model observations
     * @return A Vector which is the observations of the model
     */
    public Vector<String> getObservations() { return observations; }

    /**
     * Set the number of observations in the model
     * @param numberOfObservations An integer that specifies the number of observations in the model
     */
    public void setNumberOfObservations(int numberOfObservations) {
        this.numberOfObservations = numberOfObservations;
    }

    /**
     * Get the initial probability vector of the states
     * @return Hashtable that is the initial probability vector of the states
     */
    public Hashtable<String, Double> getInitialProbabilities() {
        return initialProbabilities;
    }

    /**
     * Set the initial probability vector of the states
     * @param initialProbabilities Hashtable that is the initial probability vector of the states
     */
    public void setInitialProbabilities(Hashtable<String, Double> initialProbabilities) {
        this.initialProbabilities = initialProbabilities;
    }

    /**
     * Get the transition matrix between the states
     * @return Hashtable that is the transition matrix between the states
     */
    public Hashtable<Pair<String, String>, Double> getTransitionMatrix() {
        return transitionMatrix;
    }

    /**
     * Set the transition matrix between the states
     * @param transitionMatrix Hashtable that is the transition matrix between the states
     */
    public void setTransitionMatrix(Hashtable<Pair<String, String>, Double> transitionMatrix) {
        this.transitionMatrix = transitionMatrix;
    }

    /**
     * Get the emission matrix between the states and the observations
     * @return Hashtable that is the emission matrix between the states and the observations
     */
    public Hashtable<Pair<String, String>, Double> getEmissionMatrix() {
        return emissionMatrix;
    }

    /**
     * Set the emission matrix between the states and the observations
     * @param emissionMatrix Hashtable that is the emission matrix between the states and the observations
     */
    public void setEmissionMatrix(Hashtable<Pair<String, String>, Double> emissionMatrix) {
        this.emissionMatrix = emissionMatrix;
    }

    /**
     *
     * @param firstState A string that is a state in the model
     * @param secondState A string that is a state in the model
     * @return A Double that is the transition value between the 2 states
     */
    public Double getTransitionValue(String firstState, String secondState) {
        return this.transitionMatrix.get(new Pair<String, String>(firstState, secondState));
    }

    /**
     *
     * @param state A string that is a state in the model
     * @param observation A string that is an observation in the model
     * @return A Double that is the value of the emission between the state and the observation
     */
    public Double getEmissionValue(String state, String observation) {
    	return this.emissionMatrix.get(new Pair<String, String>(state, observation));
    }

    /**
     *
     * @param state A string that is a state in the model
     * @return A Double that is the initial probability value of the state
     */
    public Double getInitialProbability(String state) {
        return this.initialProbabilities.get(state);
    }
    

    /**
     * Get the most optimal path for states to emit the given observations
     * @param states A Vector which is the model states
     * @param observations A Vector which represents the observations
     * @return A String which holds the optimal path and the total cost
     */
    public String getOptimalStateSequenceUsingViterbiAlgorithm(Vector<String> states, Vector<String> observations) {
        String path = "";
        Vector<Hashtable<String, Double>> dpTable = new Vector<Hashtable<String, Double>>();
        Hashtable<String, Double> statesProbabilities = new Hashtable<String, Double>();
        Hashtable<String, Double> priorProbabilities = new Hashtable<String, Double>();

        for (int i = 0; i < observations.size(); i++) {
            if (i == 0) {
                for (String state : states) {
                    double initialProbability = this.getInitialProbability(state);
                    Double emissionProbabilityObject = this.getEmissionValue(state, observations.get(i));
                    if(emissionProbabilityObject!=null){
                        double emissionProbability = getEmissionValue(state,observations.get(i));
                        statesProbabilities.put(state, Math.log(initialProbability) + Math.log(emissionProbability));
                    }
                    else{
                        statesProbabilities.put(state, 0.0);     
                    }

                }
            } else {
                for (String state : states) {
                    Double emissionProbabilityObject = this.getEmissionValue(state, observations.get(i));
                    if(emissionProbabilityObject!=null){
                    double bestProbability = -100000;
                    double emissionProbability = emissionProbabilityObject.doubleValue();
                    for (String prevState : priorProbabilities.keySet()) {
                        double transitionProbability = this.getTransitionValue(prevState, state);
                        double accumulate = priorProbabilities.get(prevState) + Math.log(emissionProbability) + Math.log(transitionProbability);

                        if (accumulate > bestProbability)
                            bestProbability = accumulate;
                    }
                    statesProbabilities.put(state, bestProbability);
                    } //END if emission probability not null
                    else{
                        statesProbabilities.put(state, 0.0);                   	
                    }
                }
            }

            dpTable.add((Hashtable<String, Double>)statesProbabilities.clone());
            priorProbabilities = (Hashtable<String, Double>) statesProbabilities.clone();
        }

        Hashtable<String, Double> lastColumn = dpTable.get(dpTable.size() - 1);
        double totalCost = -1000000;

        for (String item : lastColumn.keySet()) {
            if (lastColumn.get(item) > totalCost) {
                totalCost = lastColumn.get(item);
            }
        }

        for (Hashtable<String, Double> column : dpTable) {
            double costPerColumn = -1000000;
            String targetState = "";
            for (String state : column.keySet()) {
                if (column.get(state) > costPerColumn) {
                    costPerColumn = column.get(state);
                    targetState = state;
                }
            }
            path += targetState + " -> ";
        }

        path += "END with total cost = " + totalCost;

        return path;
    }
}