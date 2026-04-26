#include <complex>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "SimTKcommon/Testing.h"

#include "SimTKcommon.h"

using namespace SimTK;

//------------------------------------------------------------------------------
// TEST UTILITIES & MOCK TYPES
//------------------------------------------------------------------------------

template <class T>
class OtherArray_ : public Array_<T> {
    public:
    using size_type = typename Array_<T>::size_type;

    OtherArray_()
        : Array_<T>() {
    }
    OtherArray_(size_type n, const T& v)
        : Array_<T>(n, v) {
    }
};

SimTK_DEFINE_UNIQUE_INDEX_TYPE(TestIx);

/**
 * This index type has a max size of 4 for testing out-of-space checks.
 * It is used here to verify custom NiceTypeName specializations.
 */
class SmallIx {
    public:
    SmallIx()
        : ix(static_cast<unsigned char>(0xff)) {
    }
    explicit SmallIx(unsigned char i)
        : ix(i) {
    }

    auto operator++() -> SmallIx& {
        assert(ix < max_size());
        ++ix;
        return *this;
    }

    auto operator++(int) -> SmallIx {
        assert(ix < max_size());
        const SmallIx x = *this;
        ++ix;
        return x;
    }

    auto operator--() -> SmallIx& {
        assert(ix > 0);
        --ix;
        return *this;
    }

    auto operator--(int) -> SmallIx {
        assert(ix > 0);
        const SmallIx x = *this;
        --ix;
        return x;
    }

    explicit operator unsigned char() const {
        return ix;
    }
    using size_type = unsigned char;
    using difference_type = signed char;
    static auto max_size() -> size_type {
        return 4;
    }

    private:
    unsigned char ix;
};

/**
 * This class tracks constructor and destructor calls to ensure
 * types used in containers are handled correctly.
 */
class Counter {
    public:
    Counter()
        : count(0) {
    }
    auto operator=(int i) -> Counter& {
        count = i;
        return *this;
    }
    auto operator++() -> Counter& {
        ++count;
        return *this;
    }
    auto operator++(int) -> Counter {
        const Counter c = *this;
        ++count;
        return c;
    }
    auto reset() -> Counter& {
        count = 0;
        return *this;
    }
    explicit operator int() const {
        return count;
    }

    private:
    mutable int count;
};

inline auto operator<<(std::ostream& o, const Counter& c) -> std::ostream& {
    return o << static_cast<int>(c);
}

template <class T>
struct Count {
    Count() {
        ++defCtor;
    }
    Count(const Count& c)
        : val(c.val) {
        ++copyCtor;
    }
    auto operator=(const Count& c) -> Count& {
        val = c.val;
        ++copyAssign;
        return *this;
    }
    ~Count() {
        ++dtor;
    }

    Count(const T& t)
        : val(t) {
        ++initCtor;
    }
    auto operator=(const T& t) -> Count& {
        val = t;
        ++initAssign;
        return *this;
    }

    enum Color {
        Red,
        Green,
        Blue
    };
    enum class Letter {
        A,
        B,
        C
    };

    auto operator==(const Count& other) const -> bool {
        return val == other.val;
    }
    auto operator!=(const Count& other) const -> bool {
        return val != other.val;
    }

    T val;

    static auto reset() -> void {
        defCtor = 0;
        initCtor = 0;
        copyCtor = 0;
        initAssign = 0;
        copyAssign = 0;
        dtor = 0;
    }

    static Counter defCtor;
    static Counter initCtor;
    static Counter copyCtor;
    static Counter initAssign;
    static Counter copyAssign;
    static Counter dtor;
};

template <class T>
Counter Count<T>::defCtor;
template <class T>
Counter Count<T>::initCtor;
template <class T>
Counter Count<T>::copyCtor;
template <class T>
Counter Count<T>::initAssign;
template <class T>
Counter Count<T>::copyAssign;
template <class T>
Counter Count<T>::dtor;

namespace LocalNS {
enum MyEnum {
    One,
    Two
};
enum class YourEnumClass {
    Three,
    Four
};
} // namespace LocalNS

//------------------------------------------------------------------------------
// SPECIALIZATIONS
//------------------------------------------------------------------------------

namespace SimTK {
template <>
struct NiceTypeName<SmallIx> {
    static auto name() -> const char* {
        return "CustomSmallIxName";
    }
    static auto namestr() -> const std::string& {
        static const std::string ns(name());
        return ns;
    }
    static auto xmlstr() -> const std::string& {
        return namestr();
    }
};

template <class T>
struct NiceTypeName<OtherArray_<T>> {
    static auto name() -> const char* {
        return typeid(OtherArray_<T>).name();
    }
    static auto namestr() -> const std::string& {
        static const std::string ns("OtherArray_<" + NiceTypeName<T>::namestr() + ">");
        return ns;
    }
    static auto xmlstr() -> const std::string& {
        static const std::string xs = encodeTypeNameForXML(namestr());
        return xs;
    }
};
} // namespace SimTK

//------------------------------------------------------------------------------
// TEST CASES
//------------------------------------------------------------------------------

/**
 * Standardize "unsigned int" to "unsigned"; get rid of extra spaces.
 * OSX's stl likes to throw in these extra namespaces.
 */
TEST(SimTKCommon_NiceTypeName, Canonicalize_StandardizesMangledStrings) {
    EXPECT_EQ(canonicalizeTypeName("class std :: vector < unsigned int >"), "std::vector<unsigned>");
    EXPECT_EQ(canonicalizeTypeName("std:: __1 :: __23 :: set<T>"), "std::set<T>");
    EXPECT_EQ(canonicalizeTypeName("lunch bucket"), "lunch bucket");
    EXPECT_EQ(canonicalizeTypeName("std::my__1::__23x::resigned char"), "std::my__1::__23x::resigned char");
}

TEST(SimTKCommon_NiceTypeName, Builtins_ReturnsHumanReadableStrings) {
    EXPECT_EQ(NiceTypeName<bool>::namestr(), "bool");
    EXPECT_EQ(NiceTypeName<signed char>::namestr(), "signed char");
    EXPECT_EQ(NiceTypeName<int>::namestr(), "int");
    EXPECT_EQ(NiceTypeName<unsigned int>::namestr(), "unsigned");
    EXPECT_EQ(NiceTypeName<double>::namestr(), "double");
    EXPECT_EQ(NiceTypeName<float>::namestr(), "float");
    EXPECT_EQ(NiceTypeName<std::complex<double>>::namestr(), "std::complex<double>");
}

/**
 * xmlstr should replace the brackets for templatized types to be safe for tags.
 */
TEST(SimTKCommon_NiceTypeName, XML_ProvidesValidEncodedTags) {
    EXPECT_EQ(NiceTypeName<bool>::xmlstr(), "bool");
    EXPECT_EQ(NiceTypeName<unsigned int>::xmlstr(), "unsigned");
    EXPECT_EQ(NiceTypeName<std::complex<float>>::xmlstr(), "std::complex{float}");
    EXPECT_EQ(NiceTypeName<std::complex<double>>::xmlstr(), "std::complex{double}");
}

TEST(SimTKCommon_NiceTypeName, Enums_HandlesScopedAndNamespacedNames) {
    using namespace LocalNS;
    EXPECT_EQ(NiceTypeName<MyEnum>::namestr(), "LocalNS::MyEnum");
    EXPECT_EQ(NiceTypeName<YourEnumClass>::namestr(), "LocalNS::YourEnumClass");
    EXPECT_EQ(NiceTypeName<Count<double>::Color>::namestr(), "Count<double>::Color");
}

TEST(SimTKCommon_NiceTypeName, CustomTypes_SupportsUserSpecialization) {
    EXPECT_EQ(NiceTypeName<SmallIx>::namestr(), "CustomSmallIxName");
    EXPECT_EQ(NiceTypeName<TestIx>::namestr(), "TestIx");
    EXPECT_EQ(NiceTypeName<OtherArray_<int>>::namestr(), "OtherArray_<int>");
    EXPECT_EQ(NiceTypeName<OtherArray_<SmallIx>>::namestr(), "OtherArray_<CustomSmallIxName>");
}

TEST(SimTKCommon_NiceTypeName, Containers_HandlesComplexTemplates) {
    EXPECT_EQ((NiceTypeName<Array_<String, char>>::namestr()), "SimTK::Array_<SimTK::String,char>");
    EXPECT_EQ((NiceTypeName<Array_<String, char>>::xmlstr()), "SimTK::Array_{SimTK::String,char}");
    EXPECT_EQ((NiceTypeName<ArrayView_<int>>::namestr()), "SimTK::ArrayView_<int,unsigned>");
    EXPECT_EQ((NiceTypeName<Array_<Count<int>>>::xmlstr()), "SimTK::Array_{Count{int},unsigned}");
}

TEST(SimTKCommon_NiceTypeName, STL_StandardizesLibraryTypes) {
    EXPECT_EQ(NiceTypeName<std::string>::namestr(), "std::string");
    // Standard allocators should be expanded or handled per platform logic
    EXPECT_EQ((NiceTypeName<std::vector<int>>::namestr()), "std::vector<int,std::allocator<int>>");
}

namespace {
class SubTestIx : public TestIx {
    public:
    explicit SubTestIx(int ix)
        : TestIx(ix) {
    }
};
} // namespace

/**
 * Anonymous and function-local types have platform-dependent mangling.
 * We verify they at least return a non-empty string.
 */
TEST(SimTKCommon_NiceTypeName, Anonymous_CapturesLocalTypes) {
    struct MyLocalType {
        int i;
    };

    std::string anonName = NiceTypeName<SubTestIx>::namestr();
    std::string localName = NiceTypeName<MyLocalType>::namestr();

    EXPECT_FALSE(anonName.empty());
    EXPECT_FALSE(localName.empty());
}

auto main(int argc, char** argv) -> int {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}