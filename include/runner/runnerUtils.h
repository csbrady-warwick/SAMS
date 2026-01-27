#ifndef SAMS_RUNNERUTILS_H
#define SAMS_RUNNERUTILS_H

#include <type_traits>
#include <tuple>
#include "pp/callableTraits.h"

namespace SAMS{
    /**
     * Helper class to find out if a class has a member variable called dataPack. If so get it's type. If not return void.
     */
    /**
     * Void class when T does not have a dataPack member
     */
    template<typename T, typename = void>
    struct getDataPackType {
        using type = void;
    };
    /**
     * Specialisation when T has a dataPack member
     */
    template<typename T>
    struct getDataPackType<T, std::void_t<typename T::dataPack>> {
        using type = typename T::dataPack;
    };

    /**
     * Alias for getDataPackType::type
     */
    template<typename T>
    using getDataPackType_t = typename getDataPackType<T>::type;

    template<typename... T>
    struct uniqueTuple;

    /**
     * Base case for uniqueTuple = empty tuple
     */
    template<>
    struct uniqueTuple<> {
        using type = std::tuple<>;
    };


    /**
     * Recursive case for uniqueTuple, NOTE this will select the last occurrence of each type
     * so you cannot rely on the order of types being preserved. Also, note that void types are ignored.
     */
    template<typename T, typename... Rest>
    struct uniqueTuple<T, Rest...> {
    private:
        using restUnique = typename uniqueTuple<Rest...>::type;
        static constexpr bool isVoid = std::is_void_v<T>;
        static constexpr bool appearsLater = (std::is_same_v<T, Rest> || ...);
    public:
        using type = std::conditional_t<
            isVoid ||appearsLater,
            restUnique,
            decltype(std::tuple_cat(std::declval<std::tuple<T>>(), std::declval<restUnique>()))
        >;
    };

    /**
     * Specialisation to unpack a tuple and process its types into the uniqueTuple
     */
    template<typename... T, typename... Rest>
    struct uniqueTuple<std::tuple<T...>, Rest...> {
        using type = uniqueTuple<T..., Rest...>::type;
    };

    /**
     * Alias for uniqueTuple::type
     */
    template<typename... T>
    using uniqueTuple_t = typename uniqueTuple<T...>::type;

    /**
     * Helper class to get the datapacks from a list of simulations and create a unique tuple of datapack types
     */
    template<typename... Sims>
    struct dataPackTuple {
        public:
        using type = uniqueTuple_t<getDataPackType_t<Sims>...>;
    };

    template<typename... Sims>
    using dataPackTuple_t = typename dataPackTuple<Sims...>::type;

    /**
     * Helper class to return a tuple type lacking the first N elements
     */
    template<std::size_t N, typename T_tuple, typename=void>
    struct tupleTail;

    /**
     * Specialisation for N=0, return the full tuple
     */
    template<typename... Ts>
    struct tupleTail<0, std::tuple<Ts...>,void> {
        using type = std::tuple<Ts...>;
    };

    /**
     * Recursive case to remove the first element
     */
    template<std::size_t N, typename T_first, typename... T_rest>
    struct tupleTail<N, std::tuple<T_first, T_rest...>, std::enable_if_t<(N>0), void>> {
        static_assert(N <= sizeof...(T_rest)+1, "Error: N is larger than the tuple size in tupleTail");
        using type = typename tupleTail<N-1, std::tuple<T_rest...>>::type;
    };


    template<std::size_t N, typename T_tuple>
    using tupleTail_t = typename tupleTail<N, T_tuple>::type;

    /**
     * Fill a destination tuple selecting elements from a source by type
     */
    template<int level=0, typename T_dest, typename T_src>
    auto getFunctionParametersCore(T_src& src)
    {
        using srcType = std::decay_t<std::tuple_element_t<level, T_dest>>;
        if constexpr (level < std::tuple_size_v<T_dest>-1)
        {
            return std::tuple_cat(std::ref(std::get<srcType>(src)),
                                  getFunctionParametersCore<level+1, T_dest, T_src>(src));
        } else {
            return std::make_tuple(std::ref(std::get<srcType>(src)));
        }
    }

    /**
     * Get a tuple of function parameters from a function type and a source tuple
     */
    template<typename T_func, typename T_tuple>
    auto getFunctionParameters(T_tuple& tuple)
    {
        if constexpr (std::tuple_size_v<T_tuple> == 0){
            return std::tuple<>();
        } else {
            using T_call = far::callableTraits<T_func>::params;
            #pragma warning WIBBLE
            //return getFunctionParametersCore<0, T_call>(tuple);
            return std::tuple<>();
        }
    }

    /**
     * Get a tuple of elements from a source tuple by type
     */
    template<typename T_dest, typename T_src>
    auto getCallTupleElements(T_src& src)
    {
        if constexpr (std::tuple_size_v<T_dest> == 0 || std::tuple_size_v<T_src> == 0){
            return std::tuple<>();
        } else {
            return getFunctionParametersCore<0, T_dest, T_src>(src);
        }
    }

    /**
     * Get the index of a type in a tuple
     * If the type is not found, a static_assert false is triggered.
     */
    template<typename T, typename T_tuple, std::size_t level=0>
    struct tupleTypeIndex {
        static constexpr int64_t value = []() {
            if constexpr (level >= std::tuple_size_v<T_tuple>) {
                return -1; // Not found
            } else if constexpr (std::is_same_v<T, std::tuple_element_t<level, T_tuple>>) {
                return level;
            } else {
                return tupleTypeIndex<T, T_tuple, level + 1>::value;
            }
        }();
    };

    template<typename T, typename T_tuple>
    static constexpr std::size_t tupleTypeIndex_v = tupleTypeIndex<T, T_tuple>::value;

    //Get a type the consists of the union of two tuples
    template<typename ...Types>
    struct tupleUnion;

    template<typename... T1s, typename... T2s, typename... Rest>
    struct tupleUnion<std::tuple<T1s...>, std::tuple<T2s...>, Rest...> {
        using type = typename tupleUnion<uniqueTuple_t<std::tuple<T1s..., T2s...>>, Rest...>::type;
    };

    template<typename... T1s, typename... T2s>
    struct tupleUnion<std::tuple<T1s...>, std::tuple<T2s...>> {
        using type = uniqueTuple_t<std::tuple<T1s..., T2s...>>;
    };

    template<typename... T>
    using tupleUnion_t = typename tupleUnion<T...>::type;

} // namespace SAMS

#endif // SAMS_RUNNERUTILS_H