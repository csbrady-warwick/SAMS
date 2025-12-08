/*    Copyright 2025 SAMS Team
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifndef SAMS_VARIABLEREGISTRY_H
#define SAMS_VARIABLEREGISTRY_H

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>
#include "variableDef.h"
#include "mpiManager.h"

namespace SAMS{

    class variableRegistry{
        friend variableRegistry& getvariableRegistry();
        private:
        template<int i>
        friend struct MPIManager;
        std::unordered_map<std::string, variableDef> variables;
        std::vector<std::function<void(std::string)>> allocateCallbacks;
        std::unordered_map<std::string, variableDef>& getVariableMap() {
            return variables;
        }

        SAMS::MPIManager<MPI_DECOMPOSITION_RANK> &mpiMgr;
        SAMS::axisRegistry &axisReg;
        SAMS::memoryRegistry &memReg;

        public:

        variableRegistry(SAMS::MPIManager<MPI_DECOMPOSITION_RANK> &mpiMgr, SAMS::axisRegistry &axisReg, SAMS::memoryRegistry &memReg) : mpiMgr(mpiMgr), axisReg(axisReg), memReg(memReg) {}

        /**
         * Register a variable definition with a given name. If the variable already exists, make the definitions consistent.
         * @param name The name of the variable (Must be from list of physically meaningful names)
         * @param varDef The variable definition to register
         */
        void registerVariable(const std::string& name, const variableDef& varDef){
            auto it = variables.find(name);
            if(it == variables.end()){
                //New variable, insert
                variables.emplace(name, varDef);
                it = variables.find(name);
            }
            //Make consistent does some additional work, so must be called even for new variables
            it->second.makeConsistent(varDef);
        }

        /**
         * Syntatic sugar for registering a variable definition with a given name and dimensions. If the variable already exists, make the definitions consistent.
         * @param name The name of the variable (Must be from list of physically meaningful names)
         * @param varType The type of the variable (typeHandle)
         * @param memSpace The memory space of the variable (portableWrapper::arrayTags)
         * @param args The dimensions of the variable (dimension...)
         */
        template<typename... Args>
        void registerVariable(const std::string& name, typeHandle varType, portableWrapper::arrayTags memSpace, Args... args){
            static_assert(sizeof...(args) <= MAX_RANK, "Error: variableDef rank exceeds MAX_RANK");
            variableDef varDef(varType, memSpace, args...);
            registerVariable(name, varDef);
        }

        /**
         * Register a variable definition with a given name, but specifying the type as a template parameter.
         * @param name The name of the variable (Must be from list of physically meaningful names)
         * @param memSpace The memory space of the variable (portableWrapper::arrayTags)
         * @param args The dimensions of the variable (dimension...)
         * @tparam T The C++ type of the variable
         */
        template<typename T, typename... Args>
        void registerVariable(const std::string& name, portableWrapper::arrayTags memSpace, Args... args){
            static_assert(sizeof...(args) <= MAX_RANK, "Error: variableDef rank exceeds MAX_RANK");
            typeHandle varType = gettypeRegistry().getTypeID<T>();
            variableDef varDef(mpiMgr, axisReg, memReg, varType, memSpace, args...);
            registerVariable(name, varDef);
        }

        /**
         * Get a variable definition by name. Throws an error if the variable does not exist.
         * @param name The name of the variable
         * @return The variable definition
         */
        const variableDef& getVariable(const std::string& name) const {
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            return it->second;
        }

        /**
         * Get a variable definition by name. Throws an error if the variable does not exist.
         * @param name The name of the variable
         * @return The variable definition
         */
        variableDef& getVariable(const std::string& name) {
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            return it->second;
        }

        /**
         * Allocate memory for a specific variable
         * @param name The name of the variable
         */
        void allocate(const std::string& name){
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            SAMS::debug2 << "Allocating variable: " << name << std::endl;
            it->second.allocate();
            for(auto& callback : allocateCallbacks){
                callback(name);
            }
        }

        /**
         * Allocate memory for all registered variables
         */
        void allocateAll(){
            for(auto& [name, var] : variables){
                SAMS::debugAll3 << "Setting up allocation for variable: " << name << std::endl;
                allocate(name);
            }
        }


        /**
         * Deallocate memory for a specific variable
         * @param name The name of the variable
         */
        void deallocate(const std::string& name){
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            it->second.deallocate();
        }

        /**
         * Deallocate memory for all registered variables
         */
        void deallocateAll(){
            for(auto& [name, var] : variables){
                deallocate(name);
            }
        }

        /**
         * Finalize all variables in the registry
         */
        void finalize(){
            deallocateAll();
            variables.clear();
        }

        /**
         * Add a callback function to be called whenever a variable is allocated
         * @param callback The callback function taking the variable name as a parameter
         */
        void addAllocateCallback(std::function<void(std::string)> callback){
            allocateCallbacks.push_back(callback);
        }

        /**
         * Fill an internal library performance portable array from the description of a variable in the registry
         */
        template<typename T, int Arank , portableWrapper::arrayTags tag>
        void fillPPArray(const std::string name, portableWrapper::portableArray<T, Arank, tag> & array) const {
            const auto & varDef = getVariable(name);
            varDef.fillPPArray(array);
        }

        /**
         * Return an internal library performance portable array from the description of a variable in the registry
         */
        template<typename T, int Arank , portableWrapper::arrayTags tag>
        portableWrapper::portableArray<T, Arank, tag> getPPArray(const std::string name) const {
            const auto & varDef = getVariable(name);
            portableWrapper::portableArray<T, Arank, tag> array;
            varDef.fillPPArray(array);
            return array;
        }

        #ifdef USE_KOKKOS
        /*template<typename T, int Arank, portableWrapper::arrayTags tag>
        auto getKokkosView(const std::string name) const {
            const auto & varDef = getVariable(name);
            auto ppArray = varDef.getPPArray<T, Arank, tag>();
            return portableWrapper::kokkos::toView(ppArray);
        }*/
        #endif

        /** 
         * Do a halo exchange for a named variable
         */
        void haloExchange(const std::string &name){
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            it->second.haloExchange();
        }

        /** 
         * Do a halo exchange for a named variable on a specific axis
         */
        void haloExchange(const std::string &name, int axis){
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            it->second.haloExchange(axis);
        }

        /** 
         * Do a halo exchange for a named variable on a specific named axis
         * @param varName The name of the variable
         * @param axisName The name of the axis
         */
        void haloExchange(const std::string &varName, const std::string &axisName){
            auto it = variables.find(varName);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + varName + " not found in registry\n");
            }
            it->second.haloExchange(axisName);
        }

        /**
         * Do a halo exchange for a named variable on a specific axis and edge
         * @param name The name of the variable
         * @param axis The axis index
         * @param edgeType The edge type (lower, upper, both)
         */
        void haloExchange(const std::string &name, int axis, SAMS::domain::edges edgeType){
            auto it = variables.find(name);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + name + " not found in registry\n");
            }
            it->second.haloExchange(axis, edgeType);
        }

        /**
         * Do a halo exchange for a named variable on a specific named axis and edge
         * @param Name The name of the variable
         * @param axisName The name of the axis
         * @param edgeType The edge type (lower, upper, both)
         */
        void haloExchange(const std::string &varName, const std::string &axisName, SAMS::domain::edges edgeType){
            auto it = variables.find(varName);
            if(it == variables.end()){
                throw std::runtime_error("Error: variable " + varName + " not found in registry\n");
            }
            it->second.haloExchange(axisName, edgeType);
        }
    };

};

#endif //SAMS_VARIABLEREGISTRY_H