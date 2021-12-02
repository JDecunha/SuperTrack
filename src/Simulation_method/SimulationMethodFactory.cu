#include "SimulationMethodFactory.hh"
#include <iostream>

void SimulationMethodFactory::AddSimulationMethod(const std::string& name, SimulationMethod* (*constructorPointer)(const INIReader&))
{
	_nameConstructorPairs.push_back(std::make_pair(name,constructorPointer));
}

SimulationMethod* SimulationMethodFactory::Construct(const INIReader& reader)
{
	//get the simulationmethod name from the .ini file
	std::string simulationMethodName = reader.Get("Simulation","Method","");

	for (auto simMethod : _nameConstructorPairs)
	{
		if (simMethod.first == simulationMethodName)
		{
			return simMethod.second(reader);
		}
		else
		{
			std::cout << "No simulation method matching name: " << simulationMethodName << " found. Program aborting." << std::endl;
			abort();
		}
	}
}