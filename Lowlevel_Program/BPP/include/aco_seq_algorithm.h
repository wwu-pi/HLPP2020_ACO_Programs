// ===  FUNCTION HEaders ======================================================================
//      file: aco_seq_algorithm
//      author: Breno Menezes
//      Function headers
// =====================================================================================

	double run_aco_bpp(int n_ant, int n_iterations, int problem_id, int isPalma);

	double init (int nAnts, int problem, int iterations, bool is_palma);

	void initializePheromoneMatrix(int n_objects, double* phero);

	void packing();

	void pheromone_evaporation();

	void pheromone_deposit();

	void readBPPFileProperties(int problem, double* n_objects_type, double* bin_capacity);

	void readBPPFile(int problem, double* n_objects_type, double* n_objects_total, double* bin_capacity, double* items);
	
	//Print functions

	void printPOSITIONS();

	void printPHEROMONES ();

	void printGRAPH ();

	void printRESULTS ();
