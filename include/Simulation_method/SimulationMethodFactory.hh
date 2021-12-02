#pragma once

#include "SimulationMethod.hh"
#include "INIReader.h"
#include <vector>
#include <utility> //for pair, and std::make_pair
#include <string>

//Guide I used to develop the singleton
//https://stackoverflow.com/questions/1008019/c-singleton-design-pattern

class SimulationMethodFactory
{
	public:

		//Proper method to instantiate singleton
		static SimulationMethodFactory& GetInstance()
		{
			static SimulationMethodFactory singleton;

			return singleton;
		}

		//Delete the copy and = constructors
		SimulationMethodFactory(SimulationMethodFactory const&) = delete;
		void operator=(SimulationMethodFactory const&) = delete;

		//Function to add a new simulation method to the factory, just give the name and a pointer to an appropriate constructor
		void AddSimulationMethod(const std::string& name, SimulationMethod* (*constructorPointer)(const INIReader&));
		SimulationMethod* Construct(const INIReader& reader);

	private:
		SimulationMethodFactory() {} //private constructor, defined inline

		//A vector of pairs, containing the simulation method name and a constructor pointer
		std::vector<std::pair<std::string,SimulationMethod* (*)(const INIReader&)>> _nameConstructorPairs;
};