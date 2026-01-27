#ifndef SAMS_RUNNER_H
#define SAMS_RUNNER_H

#include "runnerUtils.h"
#include "pp/callableTraits.h"
#include <tuple>
#include <sstream>

namespace SAMS
{

    #define HAS_X(X) \
        private: \
        template<int level> \
        struct hasMethod_##X { \
            using type = std::tuple_element_t<level, T_combined>; \
            using params = far::callableTraits<decltype(&type::X)>::params; \
            static constexpr bool value = far::callableTraits<decltype(&type::X)>::value; \
        };

    #define HAS_X_TEMPLATE(X) \
        private: \
        template<int level, typename... T> \
        struct hasMethodTemplate_##X { \
            using type = std::tuple_element_t<level, T_combined>; \
            using params = far::callableTraits<decltype(&type::template X<T...>)>::params; \
            static constexpr bool value = far::callableTraits<decltype(&type::template X<T...>)>::value; \
        };

     /**
     * Macro to get a type and value for a named parameter which may or may not be defined
     * If it is not defined then the type is void and the value is an integer with value 0
     */
    #define HAS_X_PARAM(X)\
    private: \
        template<typename T, typename = void>\
        struct hasParam_##X {\
            using type = void;\
            static constexpr int value = 0;\
        };\
        template<typename T>\
        struct hasParam_##X<T, std::void_t<decltype(T::X)>> {\
            using type = decltype(T::X);\
            static constexpr type value = T::X;\
        }; \
        template<bool ifActive = true, bool recursive = true, int level=0>\
        T_sizeType countParam_##X(){\
            T_sizeType count;\
            if constexpr(!std::is_void_v<typename hasParamType_##X<std::tuple_element_t<level, T_combined>>::type>){\
                using type = typename hasParamType_##X<std::tuple_element_t<level, T_combined>>::type;\
                static_assert(std::is_convertible_v<type, REQ>, "Error: Parameter " #X " does not match the required type."); \
                T_sizeType shouldCount = ifActive ? simulationActiveFlags[level] : true;\
                count = shouldCount ? 1:0;\
            }\
            if constexpr (recursive && level < sizeof...(Sims)-1){\
                return count + countParam_##X<ifActive, recursive, level+1>();\
            } else {\
                return count;\
            }\
        }


    /**
     * Macro to get a type with a specified type requirement and value for a named parameter which may or may not be defined
     * If it is not defined then the type is void and the value is an integer with value 0
     * The type requirement is enforced using a static_assert
     */
    #define HAS_X_PARAM_WITH_TYPE(X, REQ)\
    private: \
        template<typename T, typename = void>\
        struct hasParamType_##X {\
            using type = void;\
            static constexpr int value = 0;\
        };\
        template<typename T>\
        struct hasParamType_##X<T, std::void_t<decltype(T::X)>> {\
            using type = decltype(T::X);\
            static constexpr type value = T::X;\
        }; \
        template<bool ifActive = true, bool recursive = true, int level=0>\
        T_sizeType countParam_##X(){\
            T_sizeType count;\
            if constexpr(!std::is_void_v<typename hasParamType_##X<std::tuple_element_t<level, T_combined>>::type>){\
                using type = typename hasParamType_##X<std::tuple_element_t<level, T_combined>>::type;\
                static_assert(std::is_convertible_v<type, REQ>, "Error: Parameter " #X " does not match the required type."); \
                T_sizeType shouldCount = ifActive ? simulationActiveFlags[level] : true;\
                count = shouldCount ? 1:0;\
            }\
            if constexpr (recursive && level < sizeof...(Sims)-1){\
                return count + countParam_##X<ifActive, recursive, level+1>();\
            } else {\
                return count;\
            }\
        }

    //Check if a simulation has a method X, and if so call it with arguments taken from dataPacks
    #define CALL_X(X)\
        HAS_X(X) \
        public: \
        template<bool ifActive = true, bool recursive = true, int level=0>\
        auto callCore_##X(){\
            if constexpr(hasMethod_##X<level>::value){\
                using type = std::tuple_element_t<level, T_combined>;\
                using fType = decltype(&type::X);\
                auto cTuple = getFunctionParameters<fType>(runnerData);\
            }\
            if constexpr (recursive && level < sizeof...(Sims)-1){\
                callCore_##X<ifActive, recursive, level+1>();\
            }\
        }\
        template<bool ifActive = true, bool recursive = true, int level=0>\
        T_sizeType getCount_##X(){\
            T_sizeType count;\
            if constexpr(hasMethod_##X<level>::value){\
                using type = std::tuple_element_t<level, T_combined>;\
                using fType = decltype(&type::X);\
                T_sizeType shouldCount = ifActive ? simulationActiveFlags[level] : true;\
                count = shouldCount ? 1:0;\
            }\
            if constexpr (recursive && level < sizeof...(Sims)-1){\
                return count + getCount_##X<ifActive, recursive, level+1>();\
            } else {\
                return count;\
            }\
        }

               /* 
               auto cTuple = getFunctionParameters<fType>(runnerData);\
               bool shouldRun = ifActive ? simulationActiveFlags[level] : true;\
                if (shouldRun) {\
                    return std::apply([&](auto&&... args) {\
                        std::get<level>(runnerData).X(std::forward<decltype(args)>(args)...);\
                    }, cTuple);\
                } \*/


    #define CALL_X_TEMPLATE(X)\
        HAS_X_TEMPLATE(X) \
        public: \
        template<typename... T, bool ifActive = true, bool recursive = true, int level=0>\
        auto callCoreTemplate_##X(T&&... args){\
            if constexpr(hasMethodTemplate_##X<level, T...>::value){\
                using rtype = hasMethodTemplate_##X<level, T...>::params; \
                using outerParams = tupleTail_t<sizeof...(T), rtype>; \
                using innerParams = decltype(std::make_tuple(std::forward<T>(args)...)); \
                bool shouldRun = ifActive ? simulationActiveFlags[level] : true;\
                if (shouldRun) {\
                    auto cTuple = std::tuple_cat( \
                        std::make_tuple(std::forward<T>(args)...), \
                        getCallTupleElements<outerParams>(runnerData) \
                    ); \
                    return std::apply([&](auto&&... cargs) {\
                        std::get<level>(runnerData).template X<T...>(std::forward<decltype(cargs)>(cargs)...);\
                    }, cTuple);\
                } \
            }\
            if constexpr (recursive && level < sizeof...(Sims)-1){\
                callCoreTemplate_##X<T... , ifActive, recursive, level+1>();\
            }\
        }

    #define FULL_CALL_X(X) \
        CALL_X(X) \
        void X(){\
            callCore_##X<>();\
        }

    #define FULL_CALL_X_TEMPLATE(X) \
        CALL_X_TEMPLATE(X) \
        template<typename... T> \
        void X(T&&... args){\
            callCoreTemplate_##X<T...>(std::forward<T>(args)...);\
        }



    /**
     * Struct that counts the number of packages that implement specific methods
     * If you add a method, you should add it here and modify setCounts accordingly
     * @note cannot be part of runner or templated since it has to be a fixed type
     */
    struct packageInfo{
        T_sizeType registerAxes;

    };

    
    /**
     * Runner class that runs a set of simulations
     */
    template<typename... Sims>
    class runner
    {
        public:
        /** Types of requested simulations. Must come first */
        using T_sims = std::tuple<Sims...>;
        harness harnessData;
        /** Internal runner data that might be requested by simulations*/
        using T_runnerData = std::tuple<harness&, packageInfo>;
        /**Data packs requested by simulations */
        using T_dataPacks = dataPackTuple_t<Sims...>;
        using T_combined = tupleUnion_t<T_sims, T_runnerData, T_dataPacks>;
        T_combined runnerData = initializeTuple<>(harnessData);

        struct internalPackageInfo{
            int level = -1;
            bool isCoreSimulation = false;
            std::string name = "";
        };

        std::map<std::string, internalPackageInfo> simulationInfoMap;
        std::array<bool, sizeof...(Sims)> simulationActiveFlags; //Active isn't in the map for performance reasons

        CALL_X(initialize);//Initialize simulations
        FULL_CALL_X(registerAxes);//Register axes

        HAS_X_PARAM_WITH_TYPE(name, std::string);
        HAS_X_PARAM_WITH_TYPE(coreSimulation, bool);

        public:

        /**
         * Populate a single simulation info struct
         */
        template<int level=0>
        internalPackageInfo buildSimulationInfo(){
            internalPackageInfo info;
            info.name = hasParamType_name<std::tuple_element_t<level, T_combined>>::value;
            info.level = level;
            if constexpr (std::is_void_v<typename hasParamType_coreSimulation<std::tuple_element_t<level, T_combined>>::type> == false){
                info.isCoreSimulation = hasParamType_coreSimulation<std::tuple_element_t<level, T_combined>>::value;
            } else {
                info.isCoreSimulation = false;
            }
            return info;
        }

        /**
         * Build the simulation info map
         */
        template<int level=0>
        void buildInfo(){
            if constexpr (std::is_void_v<typename hasParamType_name<std::tuple_element_t<level, T_combined>>::type> == false){
                std::string name = hasParamType_name<std::tuple_element_t<level, T_combined>>::value;
                simulationInfoMap[name] = buildSimulationInfo<level>();
            } else {
                std::stringstream ss;
                ss << "ERROR: Simulation at index " << level << " does not have a name parameter defined.\n";
                abort(ss.str(), false);
            }
            if constexpr (level < sizeof...(Sims)-1){
                buildInfo<level+1>();
            }
        }

        /**
         * Function to initialize the combined tuple of simulations, runner data and data packs
         */
        template<int level=0>
        static auto initializeTuple(harness& h){
           constexpr std::size_t N = std::tuple_size_v<T_combined>;
            using T_current = std::tuple_element_t<level, T_combined>;

            //Case where the tuple element is just default constructible
            if constexpr (std::is_default_constructible_v<T_current>) {
                if constexpr (level + 1 < N) {
                    return std::tuple_cat(
                        std::make_tuple(T_current{}),
                        initializeTuple<level + 1>(h)
                    );
                } else {
                    return std::make_tuple(T_current{});
                }
                //Case where the tuple element is a reference to the harness
            } else if constexpr (std::is_reference_v<T_current>) {
                if constexpr (level + 1 < N) {
                    return std::tuple_cat(
                        std::tuple<T_current>(h),
                        initializeTuple<level + 1>(h)
                    );
                } else {
                    return std::tuple<T_current>(h);
                }
                //Case where the tuple element is constructible from a harness reference
            } else if constexpr (std::is_constructible_v<T_current, harness&>) {
                if constexpr (level + 1 < N) {
                    return std::tuple_cat(
                        std::make_tuple(T_current{h}),
                        initializeTuple<level + 1>(h)
                    );
                } else {
                    return std::make_tuple(T_current{h});
                }
            } else {
                static_assert(portableWrapper::alwaysFalse<T_current>::value,
                              "Error: Simulation or data pack type is not default constructible or constructible from a harness reference.");
            }
        }

        public:

        /**
         * Abort the runner, dealing with MPI as needed
         */
        void abort(const std::string& message = "", bool localError = false){
            getHarness().abort(message, localError);
        }

        /**
         * Is T a simulation type managed by this runner
         * @tparam T The type to check
         */
        template<typename T>
        constexpr bool isSimulationType(){
            constexpr int64_t index = tupleTypeIndex_v<T, T_combined>;
            if constexpr (index < 0 || index >= sizeof...(Sims)){
                return false;
            } else {
                return true;
            }
        }

        /**
         * Is T a runner data type managed by this runner
         * @tparam T The type to check
         */
        template<typename T>
        constexpr bool isRunnerDataType(){
            constexpr int64_t index = tupleTypeIndex_v<T, T_combined>;
            if constexpr (index < 0 || index < sizeof...(Sims) || index >= sizeof...(Sims) + std::tuple_size_v<T_runnerData>){
                return false;
            } else {
                return true;;
            }
        }

        /**
         * is T a data pack type managed by this runner
         * @tparam T The type to check
         */
        template<typename T>
        constexpr bool isDataPackType(){
            constexpr int64_t index = tupleTypeIndex_v<T, T_combined>;
            if constexpr (index < 0 || index >= sizeof...(Sims) + std::tuple_size_v<T_runnerData>){
                return false;
            } else {
                return true;
            }
        }

        /**
         * Get a physics package from the runner
         * @tparam T_sim The simulation type
         * @return A reference to the simulation
         */
        template<typename T_sim>
        T_sim& getSimulation()
        {
            static_assert(isSimulationType<T_sim>(), "Error: T_sim is not a simulation type managed by this runner.");
            //Should I check that this is a simulation type, not a data pack or runner data?
            return std::get<T_sim>(runnerData);
        }

        /**
         * Get the harness from the runner
         * @return A reference to the harness
         */
        harness& getHarness()
        {
            return std::get<harness&>(runnerData);
        }

        /**
         * Get a data pack from the runner
         * @tparam T_dataPack The data pack type
         * @return A reference to the data pack
         */
        template<typename T_dataPack>
        T_dataPack& getDataPack()
        {
            static_assert(isDataPackType<T_dataPack>(), "Error: T_dataPack is not a data pack type managed by this runner.");
            return std::get<T_dataPack>(runnerData);
        }


        //Initialise the runner
        void initialize(int &argc, char** &argv)
        {
            getHarness().initialize(argc, argv);
            //Turn off all simulations
            std::fill(simulationActiveFlags.begin(), simulationActiveFlags.end(), false);
            //Build the simulation info map
            buildInfo();
        }

        // Activate a simulation by name
        void activateSimulation(const std::string& simName)
        {
            auto it = simulationInfoMap.find(simName);
            if (it != simulationInfoMap.end()){
                simulationActiveFlags[it->second.level] = true;
            } else {
                std::stringstream ss;
                ss << "ERROR: Simulation with name " << simName << " not found in runner.\n";
                abort(ss.str(), false);
            }
        }


        //Iniitialize all active simulations
        void initializeSimulations()
        {
            buildInfo();
            callCore_initialize();
        }

    };
    
} // namespace SAMS


#endif // SAMS_RUNNER_H