#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include "pico/stdlib.h"
#include "pico/rand.h"

#define CODE_LENGTH	4
#define COLUMN_LENGTH 7
#define COLUMNS 6

using namespace std;

/**----------------------------------------------
 *                     ABOUT
 * @author      : Uriel Hansen
 * @email       : urielhansen2005@gmail.com
 * @repo        : NA
 * @createdOn   : 1/23/2025
 * @description : Data generator for a keypad code.
 *---------------------------------------------**/

/**--------------------------------------------
 *              CLASS DECLARATIONS
 *---------------------------------------------**/

// A class to collect relevant numeric functions.
class Numeric {
    public:
        /**
		 * Returns a random integer between min and max.
		 * @param min Minimum number generated (inclusive).
		 * @param max Maximum number generated (exclusive).
		 * @return The randomly generated integer.
		 **/
        static uint8_t random (uint8_t, uint8_t);
        /**
		 * Returns the factorial of the input number.
		 * @param input Input number.
		 * @return Factorial of input number.
		**/
        static uint16_t factorial (uint8_t);
        /**
		 * Returns the number of configurations of choosing k items from set of n items.
		 * @param n Number of total items.
		 * @param k Number of items selected.
		 * @return Number of configurations possible.
		**/
        static uint8_t nCk (uint8_t, uint8_t);
};

/*
 * A class to choose a random category, predetermined number of elements from 
 * the category, and random permutation order to display the elements.
 * @tparam T Datatype/class of elements.
*/
template<typename T> class OrderedCategoryRandomizer {
    private:
		/**
		 * Vector of pointers to configuration data.
		 * 
		 * NOTE: Configurations are comprised of all possible ascending orderings of (codeLength) number
		 * of indices in the range 0 to catyegorySize. Configuration data is stored as an array of indices.
		**/
        vector<uint8_t*> configurations;
		/**
		 * Vector of pointers to permutation data. 
		 * 
		 * NOTE: Permutations are comprised of all possible orderings of (codeLength) 
		 * number of elements. Permutation data is stored as an array of indices.
		**/
        vector<uint8_t*> permutations;
		
		// Number of categories.
		uint8_t categoryCt;
		// Size of categories.
		uint8_t categorySize;
		// Length of code to be generated.
		uint8_t codeLength;
		vector<T*> layout;

		// Number of configurations generated.
		uint8_t configurationCt;
		// Number of permutations generated.
		uint8_t permutationCt;
		
		/**
		 * Pre-generates all of the possible configurations (ran from constructor).
		 * @return void
		**/
		void generate_configurations();
		/**
		 * Pre-generates all of the possible permutations (ran from constructor).
		 * @return void
		**/
		void generate_permutations();
		/**
		 * Iterative step in configuration generation.
		 * @param length Current length of code in recursion.
		 * @param configData Array of current elements of code of configuration.
		 * @return void
		**/
		void generate_config_iterate(uint8_t, uint8_t[]);
		/**
		 * Iterative step in permutation generation.
		 * @param length Current length of code in recursion.
		 * @param permuteData Array of current elements of code of permutation.
		 * @return void
		**/
		void generate_permut_iterate(uint8_t, uint8_t[]);
	public:
		/**
		 * Constructor
		 * @param categoryCt Number of categories provided.
		 * @param categorySize Size of the categories provided.
		 * @param codeLength Length of code to be generated.
		 * @param layout A vector consisting of pointers to arrays of elements which represent the categories.
		**/
		OrderedCategoryRandomizer(uint8_t, uint8_t, uint8_t, vector<T*>);
		
		/**
		 * Returns a configuration by numeric ID.
		 * @param ID ID number of configuration.
		 * @return Pointer to array of configuration data.
		**/
		uint8_t* configuration(uint8_t);
		/**
		 * Returns a permutation by numeric ID.
		 * @param ID ID number of permutation.
		 * @return Pointer to array of permutation data.
		**/
		uint8_t* permutation(uint8_t);
		/**
		 * Returns the code length created from the randomizer.
		 * @return Length of code to be generated.
		**/
		uint8_t code_length ();
		
		/**
		 * Returns the category count of the randomizer.
		 * @return Number of categories.
		**/
		uint8_t category_ct ();
		/**
		 * Returns the category size of the randomizer.
		 * @return Size of categories.
		**/
		uint8_t category_sz ();
		/**
		 * Returns the number of configurations in the randomizer.
		 * @return Number of configurations generated.
		**/
		uint8_t configuration_ct ();
		/**
		 * Returns the number of permutations in the randomizer.
		 * @return Number of permutations generated.
		**/
		uint8_t permutation_ct ();

		/**
		 * Returns the element of the category with the given ID at the given index.
		 * @param category ID number of category.
		 * @param index Index of element in the category.
		 * @return The element at the index position in the category with the given ID.
		**/
		T get(uint8_t, uint8_t);

		//Stores data for a particular generated code.
		struct Instance {
			//Category ID of the instance.
			uint8_t categoryID;
			//Configuration ID of the instance.
			uint8_t configID;
			//Permutation ID of the instance
			uint8_t permutationID;

			// Pointer to array of elements in display order.
			T* display;
			// Pointer to array of elements in original (correct) order.
			T* key;

			Instance(uint8_t cat, uint8_t cfg, uint8_t prm, T* dpl, T* key)
				: categoryID(cat), configID(cfg), permutationID(prm), display(dpl), key(key) {};
		};

		/**
		 * Returns a set of randomized data stored in a struct.
		 * @return Struct of randomized data of datatype OrderedCategoryRandomizer::Instance<T>.
		**/
		struct Instance random();
};

/**--------------------------------------------
 *             MAIN LOOP (testing)
 *---------------------------------------------**/

int main () {
	stdio_init_all ();

	vector<char*> layout;
	layout.push_back (new char[] { 'A', 'B', 'C', 'D', 'E', 'F', 'G' });
	layout.push_back (new char[] { 'H', 'A', 'G', 'I', 'J', 'F', 'K' });
	layout.push_back (new char[] { 'L', 'M', 'I', 'N', 'O', 'C', 'J' });
	layout.push_back (new char[] { 'P', 'Q', 'R', 'E', 'N', 'K', 'S' });
	layout.push_back (new char[] { 'T', 'S', 'R', 'U', 'Q', 'V', 'W' });
	layout.push_back (new char[] { 'P', 'H', 'X', 'Y', 'T', 'Z', '0' });

	while (true) {
		getchar();

		OrderedCategoryRandomizer<char> cfg = OrderedCategoryRandomizer (COLUMNS, COLUMN_LENGTH, CODE_LENGTH, layout);
		OrderedCategoryRandomizer<char>::Instance inst = cfg.random ();
	}

    return 0;
}

/**--------------------------------------------
 *              CLASS DEFINITIONS
 *---------------------------------------------**/

uint8_t Numeric::random (uint8_t min, uint8_t max) {
	return (uint8_t)(min + (get_rand_32() % (max - min)));
}

uint16_t Numeric::factorial (uint8_t n) {
    assert (n >= 0);

    uint16_t val = 1;
    for (uint16_t i = 1; i <= n; i++)
        val *= i;
    return val;
}

uint8_t Numeric::nCk (uint8_t n, uint8_t k) {
    assert ((n >= 0) * (k >= 0));

    if (k > n) { return 0; }
	return (uint8_t)(Numeric::factorial (n) / (Numeric::factorial (k) * Numeric::factorial (n - k)));
}

template<typename T> OrderedCategoryRandomizer<T>::OrderedCategoryRandomizer(uint8_t categoryCt, uint8_t categorySize, uint8_t codeLength, vector<T*> layout) {
	this->categoryCt = categoryCt;
	this->categorySize = categorySize;
	this->codeLength = codeLength;
	this->layout = layout;

	generate_configurations ();
	generate_permutations ();
}

template<typename T> void OrderedCategoryRandomizer<T>::generate_config_iterate (uint8_t length, uint8_t configData[]) {
	assert ((length >= 0) * (length <= codeLength));

	uint8_t i, mins;

	// End of recursion, full config is generated.
	if (length == codeLength) {
		// Creates copy of configData to push into configurations to avoid overwriting of data in further iterations.
		uint8_t *tmp = new uint8_t[codeLength];
		for (i = 0; i < length; i++)
			tmp[i] = configData[i];

		configurations.push_back(tmp);
		return;
	}

	// Minimum element is (previous element + 1), unless there is no previous element, then 0.
	mins = (length == 0 ? 0 : configData[length - 1] + 1);
	for (i = mins; i < categorySize - codeLength + length + 1; i++) {
		// Iterate
		configData[length] = i;
		generate_config_iterate (length + 1, configData);
	}
}

template<typename T> void OrderedCategoryRandomizer<T>::generate_permut_iterate (uint8_t length, uint8_t permuteData[]) {
	assert ((length >= 0) * (length <= codeLength));

	uint8_t i, j;
	bool flag;

	// End of recursion, full permutation is generated.
	if (length == codeLength) {
		// Creates copy of permuteData to push into permutations to avoid overwriting of data in further iterations.
		uint8_t *tmp = new uint8_t[codeLength];
		for (i = 0; i < length; i++)
			tmp[i] = permuteData[i];

		permutations.push_back(tmp);
		return;
	}

	for (i = 0; i < codeLength; i++) {
		// Confirms i is not included in the current permuteData.
		flag = 0;
		for (j = 0; j < length; j++)
			{ if (permuteData[j] == i) flag = 1; }
		if (flag) continue;

		// Iterate
		permuteData[length] = i;
		generate_permut_iterate (length + 1, permuteData);
	}
}

template<typename T> void OrderedCategoryRandomizer<T>::generate_configurations () {
	// Calculates number of configurations to be generated.
	this->configurationCt = Numeric::nCk (categorySize, codeLength);

	// Begins iterative loop of configuration generation.
	uint8_t data[codeLength];
	this->generate_config_iterate (0, data);
}

template<typename T> void OrderedCategoryRandomizer<T>::generate_permutations () {
	// Calculates number of permutations to be generated.
	this->permutationCt = Numeric::factorial (codeLength);

	// Begins iterative loop of permutation generation.
	uint8_t data[codeLength];
	this->generate_permut_iterate (0, data);
}

template<typename T> uint8_t* OrderedCategoryRandomizer<T>::configuration(uint8_t i) {
	assert((i >= 0) * (i < this->configurationCt));
	return this->configurations.at (i);
}

template<typename T> uint8_t* OrderedCategoryRandomizer<T>::permutation(uint8_t i) {
	assert((i >= 0) * (i < this->permutationCt));
	return this->permutations.at (i);
}

template<typename T> uint8_t OrderedCategoryRandomizer<T>::code_length () { return this->codeLength; }
template<typename T> uint8_t OrderedCategoryRandomizer<T>::configuration_ct () { return this->configurationCt; }
template<typename T> uint8_t OrderedCategoryRandomizer<T>::permutation_ct () { return this->permutationCt; }
template<typename T> uint8_t OrderedCategoryRandomizer<T>::category_ct() { return this->categoryCt; }
template<typename T> uint8_t OrderedCategoryRandomizer<T>::category_sz() { return this->categorySize; }

template<typename T> T OrderedCategoryRandomizer<T>::get(uint8_t category, uint8_t index) {
	assert ((category >= 0) * (category < this->category_ct ()) * (index >= 0) * (index < this->category_sz ()));
	return this->layout.at (category)[index];
}

template<typename T> struct OrderedCategoryRandomizer<T>::Instance OrderedCategoryRandomizer<T>::random() {
	uint8_t i;
	
	// Generates random category, configuration, and permutation.
	uint8_t catID = Numeric::random (0, this->category_ct());
	uint8_t cfgID = Numeric::random (0, this->configuration_ct());
	uint8_t prmID = Numeric::random (0, this->permutation_ct());

	// Initialize arrays (as pointers) for display and key.
	T* dpl = new T[this->codeLength];
	T* key = new T[this->codeLength];

	// Fills dpl and key with calculated values.
	for (i = 0; i < this->codeLength; i++) {
		key[i] = this->get(catID, this->configuration (cfgID)[i]);
		dpl[i] = this->get(catID, this->configuration (cfgID)[this->permutation (prmID)[i]]);
	}

	// Creates struct with previously generated data.
	struct OrderedCategoryRandomizer<T>::Instance inst(catID, cfgID, prmID, dpl, key);
	return inst;
}