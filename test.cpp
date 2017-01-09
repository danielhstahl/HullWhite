#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"


TEST_CASE("Test for_each_parallel", "[Functional]"){
    std::vector<int> testV={5, 6, 7};
    auto squareTestV=[](auto& val, const auto& index){
        return val*val;
    };
    REQUIRE(futilities::for_each_parallel(std::move(testV), squareTestV)==std::vector<int>({25, 36, 49}));
}
TEST_CASE("Test for_each_parallel iterator", "[Functional]"){
    auto squareTestV=[](auto& val){
        return val*val;
    };
    REQUIRE(futilities::for_each_parallel(5, 8, squareTestV)==std::vector<int>({25, 36, 49}));
}
TEST_CASE("Test for_each", "[Functional]"){
    std::vector<int> testV={5, 6, 7};
    auto squareTestV=[](const auto& val, const auto& index){
        return val*val;
    };
    REQUIRE(futilities::for_each(std::move(testV), squareTestV)==std::vector<int>({25, 36, 49}));
}
TEST_CASE("Test for_emplace_back", "[Functional]"){
    auto valTestV=[](const auto& val){
        return val;
    };
    REQUIRE(futilities::for_emplace_back(0.0, 1.0, 5, valTestV)==std::vector<double>({0, .25, .5, .75, 1}));
}
TEST_CASE("Test cumulative sum", "[Functional]"){
    std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val, const auto& index){
        return val;
    };
    REQUIRE(futilities::cumulative_sum(std::move(testV), valTestV)==std::vector<int>({5, 11, 18, 26, 35}));
}
TEST_CASE("Test sum", "[Functional]"){
    std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val, const auto& index){
        return val;
    };
    REQUIRE(futilities::sum(testV, valTestV)==35);
}
TEST_CASE("Test sum iterator", "[Functional]"){
    //std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val){
        return val;
    };
    REQUIRE(futilities::sum(5, 10, valTestV)==35);
}
TEST_CASE("Test sum nontrivial", "[Functional]"){
    //std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val){
        return val*val;
    };
    REQUIRE(futilities::sum(5, 10, valTestV)==255);
}
 
TEST_CASE("Test recurse", "[Functional]"){
    //std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val, const auto& index){
        return val*2;
    }; 
    REQUIRE(futilities::recurse(5, 1, valTestV)==pow(2, 5));
} 
TEST_CASE("Test vector recurse", "[Functional]"){
    //std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val, const auto& index){
        return std::vector<double>({val[0]*2, val[1]*2});
    }; 
    auto myTest=futilities::recurse(5, std::vector<double>({1, 1}), valTestV);
    REQUIRE(myTest[0]==pow(2, 5));
} 
 
/*TEST_CASE("Test compilation", "[Functional]"){
    std::vector<int> testV={5, 6, 7, 8, 9};
    auto valTestV=[](const auto& val, const auto& index){
        return val;
    };
    REQUIRE(futilities::cumulative_sum(testV, valTestV)==std::vector<int>({5, 11, 18, 26, 35}));
}*/
