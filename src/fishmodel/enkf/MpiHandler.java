package fishmodel.enkf;

import mpi.*;

public class MpiHandler {

    int rank, N;

    public MpiHandler(String[] args) {
        MPI.Init(args);
        rank = MPI.COMM_WORLD.Rank();
        N = MPI.COMM_WORLD.Size();
    }

    public int getRank() {
        return rank;
    }

    public int getN() {
        return N;
    }

    public double[][] gatherStateToRank0(double[][][] field) {
        if (rank==0) {
            double[][] allStates = new double[N][];
            double[] vector = Util.reshapeArray(field);
            allStates[0] = new double[vector.length];
            System.arraycopy(vector, 0, allStates[0], 0, vector.length);
            for (int i=1; i<N; i++) {
                MPI.COMM_WORLD.Recv(vector, 0, vector.length, MPI.DOUBLE, i, 1);
                allStates[i] = new double[vector.length];
                System.arraycopy(vector, 0, allStates[i], 0, vector.length);
                System.out.println("Received data from rank "+i);
            }
            return allStates;
        }
        else {
            double[] vector = Util.reshapeArray(field);
            MPI.COMM_WORLD.Send(vector, 0, vector.length, MPI.DOUBLE, 0, 1);
            return null;
        }
    }

}