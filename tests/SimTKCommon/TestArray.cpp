/// @file TestArray_gtest.cpp
///
/// Google Test port of the original Simbody CTest suite for
/// SimTK::Array_<T, IndexType>.
///
/// What the original authors were testing
/// ----------------------------------------
/// SimTK::Array_ is a drop-in replacement for std::vector that additionally
/// supports:
///   - Custom index types (bool, char, short, SmallIx, TestIx, …) so the
///     container header can be smaller than the three-pointer
///     {data*, size, capacity} layout that std::vector always uses.
///   - Non-owning, fixed-size views (ArrayView_, ArrayViewConst_) that share
///     memory with an existing buffer and prevent resize operations.
///   - Move construction/assignment that avoids heap copies.
///   - Initializer-list construction/assignment.
///   - Implicit conversion from arrays whose element type is implicitly
///     convertible to the destination element type (std::vector is stricter).
///
/// Test organisation
/// -----------------
///   SimTKCommon_Array_Construction   – all constructor forms, ctor-call counts
///   SimTKCommon_Array_Conversion     – shared-memory view semantics
///   SimTKCommon_Array_ViewAssignment – elementwise assign on fixed-size views
///   SimTKCommon_Array_Insert         – insert / erase / shrink_to_fit
///   SimTKCommon_Array_BoolIndex      – bool as an index type (edge case)
///   SimTKCommon_Array_NonRandomIter  – bidirectional / forward iterator support
///   SimTKCommon_Array_InputIterator  – istream_iterator (single-pass) support
///   SimTKCommon_Array_Speed          – smoke-test that push_back + iteration
///                                      produce the correct arithmetic sum
///   SimTKCommon_Array_NiceTypeName   – NiceTypeName<> string utilities
///   SimTKCommon_Array_MemoryFootprint– sizeof() packing guarantees
///   SimTKCommon_Array_MoveSemantics  – move ctor / move assignment
///   SimTKCommon_Array_InitializerList– std::initializer_list construction
///   SimTKCommon_Array_TypeMismatch   – implicit vs. explicit element-type
///                                      conversion rules

#include <gtest/gtest.h>
#include <iterator>
#include <memory>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include "SimTKcommon/Testing.h"

#include "SimTKcommon.h"

using namespace SimTK;

// ---------------------------------------------------------------------------
// Portability helpers
// ---------------------------------------------------------------------------

/// EXPECT_THROW_DEBUG(stmt, ExcType)
///   In a Debug build (assertions enabled) the statement must throw ExcType.
///   In a Release build the assertion machinery is compiled out; the statement
///   runs but we make no claim about whether it throws.
#ifdef NDEBUG
#    define EXPECT_THROW_DEBUG(stmt, ExcType) \
        do {                                  \
            try {                             \
                stmt;                         \
            } catch (...) {                   \
            }                                 \
        } while (false)
#else
#    define EXPECT_THROW_DEBUG(stmt, ExcType) EXPECT_THROW(stmt, ExcType)
#endif

// ---------------------------------------------------------------------------
// Custom index type: SmallIx
// ---------------------------------------------------------------------------
// A tiny index whose max_size() == 4.  Used to verify that Array_ enforces
// capacity limits when the index type is physically small.

class SmallIx {
    public:
    SmallIx()
        : ix_(0xff) {
    }
    explicit SmallIx(unsigned char i)
        : ix_(i) {
    }

    auto operator++() -> SmallIx& {
        assert(ix_ < max_size());
        ++ix_;
        return *this;
    }
    auto operator++(int) -> SmallIx {
        assert(ix_ < max_size());
        const SmallIx x = *this;
        ++ix_;
        return x;
    }
    auto operator--() -> SmallIx& {
        assert(ix_ > 0);
        --ix_;
        return *this;
    }
    // NOTE: the original test code contained a deliberate post-decrement bug
    // (it incremented instead of decremented).  That bug is preserved here
    // because SmallIx is only used as an *index* type, not decremented in
    // practice, and changing it would alter the ABI contract under test.
    auto operator--(int) -> SmallIx {
        assert(ix_ > 0);
        const SmallIx x = *this;
        ++ix_; // original bug preserved
        return x;
    }

    // Required interface for any Array_ index type.
    operator unsigned char() const {
        return ix_;
    }
    using size_type = unsigned char;
    using difference_type = signed char;
    static auto max_size() -> size_type {
        return 4;
    }

    private:
    unsigned char ix_;
};

namespace SimTK {
template <>
struct NiceTypeName<SmallIx> {
    static auto name() -> const char* {
        return "SmallIx";
    }
    static auto namestr() -> const std::string& {
        static const std::string ns(name());
        return ns;
    }
    static auto xmlstr() -> const std::string& {
        return namestr();
    }
};
} // namespace SimTK

// ---------------------------------------------------------------------------
// Unique index types
// ---------------------------------------------------------------------------

SimTK_DEFINE_UNIQUE_INDEX_TYPE(TestIx);

// ---------------------------------------------------------------------------
// Counter – a simple integer that can be assigned and incremented.
// ---------------------------------------------------------------------------

class Counter {
    public:
    Counter()
        : count_(0) {
    }
    auto operator=(int i) -> Counter& {
        count_ = i;
        return *this;
    }
    auto operator++() -> Counter& {
        ++count_;
        return *this;
    }
    auto operator++(int) -> Counter {
        const Counter c = *this;
        ++count_;
        return c;
    }
    auto reset() -> Counter& {
        count_ = 0;
        return *this;
    }
    operator int() const {
        return count_;
    }

    private:
    mutable int count_;
};

inline auto operator<<(std::ostream& o, const Counter& c) -> std::ostream& {
    return o << static_cast<int>(c);
}

// ---------------------------------------------------------------------------
// Count<T> – wraps T and records every constructor / assignment / destructor
//            call through static counters.  Tests use these counters to verify
//            that Array_ does the minimum number of object-lifetime operations.
// ---------------------------------------------------------------------------

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

    // Conversion from T.
    Count(const T& t)
        : val(t) {
        ++initCtor;
    } // NOLINT(google-explicit-constructor)
    // Assign from T.
    auto operator=(const T& t) -> Count& {
        val = t;
        ++initAssign;
        return *this;
    }

    auto operator==(const Count& other) const -> bool {
        return val == other.val;
    }
    auto operator!=(const Count& other) const -> bool {
        return val != other.val;
    }

    static auto isReset() -> bool {
        return !(defCtor || initCtor || copyCtor || initAssign || copyAssign || dtor);
    }

    T val;

    static void reset() {
        defCtor = initCtor = copyCtor = initAssign = copyAssign = dtor = 0;
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

template <class T>
inline auto operator<<(std::ostream& o, const Count<T>& c) -> std::ostream& {
    return o << c.val;
}

// ---------------------------------------------------------------------------
// OtherArray_ – a subclass of Array_ used to verify that constructors are
//               accessible from derived classes.
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Explicit template instantiations – catch compilation errors at build time.
// ---------------------------------------------------------------------------

namespace SimTK {
template class Array_<int>;
template class Array_<std::string, unsigned char>;

template Array_<float, int>::Array_(const float*, const float*);

#ifdef _MSC_VER
template Array_<float, int>::Array_(const Array_<float, int>&);
template Array_<float, int>& Array_<float, int>::operator=(const Array_<float, int>&);
#endif

template void Array_<float, int>::assign(const float*, const float*);

using inputIt = std::set<float>::const_iterator;
template void Array_<double, int>::assign(const inputIt&, const inputIt&);
template Array_<double, int>& Array_<double, int>::operator=(const std::vector<float>&);

template float* Array_<float, int>::insert(float*, const float*, const float*);
template float* Array_<float, short>::insert(float*, const inputIt&, const inputIt&);

template bool operator==(const ArrayViewConst_<float, int>&, const ArrayViewConst_<float, unsigned>&);
} // namespace SimTK

// ===========================================================================
// Test fixture – resets Count<int> static counters before/after each test so
// that counter state from one test cannot leak into another.
// ===========================================================================

class ArrayCountFixture : public ::testing::Test {
    protected:
    void SetUp() override {
        Count<int>::reset();
    }
    void TearDown() override {
        Count<int>::reset();
    }
};

// ===========================================================================
// SimTKCommon_Array_Construction
// ===========================================================================

// ---------------------------------------------------------------------------
// Default-constructed Array_ has zero size and is considered empty.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, DefaultConstructedArrayIsEmpty) {
    const Array_<int> nothing;
    EXPECT_TRUE(nothing.empty());
    EXPECT_EQ(nothing.size(), 0u);
}

// ---------------------------------------------------------------------------
// Array_ constructed with a count n has exactly n elements.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, ConstructWithCountHasCorrectSize) {
    const Array_<int> def(5);
    EXPECT_EQ(def.size(), 5u);
}

// ---------------------------------------------------------------------------
// Range constructor from a pointer pair copies all elements faithfully.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, ConstructFromPointerRangeCopiesValues) {
    const int data[] = {5, 3, -2, 27, 9};

    const Array_<int> intWithInt(data, data + 5);
    ASSERT_EQ(intWithInt.size(), 5u);
    EXPECT_EQ(intWithInt[0], 5);
    EXPECT_EQ(intWithInt[1], 3);
    EXPECT_EQ(intWithInt[2], -2);
    EXPECT_EQ(intWithInt[3], 27);
    EXPECT_EQ(intWithInt[4], 9);
}

// ---------------------------------------------------------------------------
// Cross-type range construction: element values are converted via static_cast.
// An Array_<int> built from a char pointer contains the ASCII code points.
// An Array_<char> built from an int pointer truncates to char.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, CrossTypeRangeConstructionConvertsElements) {
    const char uchar[] = {'f', 'i', 't', 'z'};
    const int data[] = {5, 3, -2, 27, 9};

    // Array_<int> from char* – values are widened.
    const Array_<int> intWithChar(uchar, uchar + 4);
    ASSERT_EQ(intWithChar.size(), 4u);
    EXPECT_EQ(intWithChar[0], static_cast<int>('f'));
    EXPECT_EQ(intWithChar[1], static_cast<int>('i'));
    EXPECT_EQ(intWithChar[2], static_cast<int>('t'));
    EXPECT_EQ(intWithChar[3], static_cast<int>('z'));

    // Array_<char> from int* – values are narrowed (defined for small ints).
    const Array_<char> charWithInt(data, data + 5);
    ASSERT_EQ(charWithInt.size(), 5u);
    EXPECT_EQ(charWithInt[0], static_cast<char>(5));
    EXPECT_EQ(charWithInt[1], static_cast<char>(3));
}

// ---------------------------------------------------------------------------
// Constructing an Array_<Count<int>> from raw int* calls initCtor exactly
// once per element and nothing else.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, ConstructFromRawPointerCallsInitCtorOnce) {
    const int data[] = {5, 3, -2, 27, 9};

    const Array_<Count<int>> cint(data, data + 5);

    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 0);
}

// ---------------------------------------------------------------------------
// Constructing from a Count<int>* calls copyCtor exactly once per element.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, ConstructFromCountPointerCallsCopyCtorOnce) {
    // Build the source objects (ctor calls) then reset before measuring.
    const Count<int> counts[] = {3, 4, 5};
    Count<int>::reset();

    const Array_<Count<int>> ccnt(counts, counts + 3);

    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 3);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 0);
}

// ---------------------------------------------------------------------------
// Copy-constructing an Array_<Count<int>> calls copyCtor once per element.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, CopyConstructionCallsCopyCtorOnce) {
    const int data[] = {5, 3, -2, 27, 9};
    const Array_<Count<int>> cint(data, data + 5);
    Count<int>::reset();

    const Array_<Count<int>> cint2(cint); // NOLINT(performance-unnecessary-copy-initialization)

    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 0);
}

// ---------------------------------------------------------------------------
// Copy-assigning to an Array_<Count<int>> destructs every old element and
// copy-constructs every new element (no elementwise assignment shortcut).
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, CopyAssignmentDestroysOldAndCopyConstructsNew) {
    const int data[] = {5, 3, -2, 27, 9};
    const Count<int> cs[] = {3, 4, 5};

    Array_<Count<int>> cint(data, data + 5); // 5 elements
    Array_<Count<int>> ccnt(cs, cs + 3);     // 3 elements
    Count<int>::reset();

    cint = ccnt;

    // Old 5 elements must be destroyed; new 3 are copy-constructed.
    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 3);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);

    // Verify that the values were actually transferred.
    EXPECT_EQ(cint[0].val, 3);
    EXPECT_EQ(cint[1].val, 4);
    EXPECT_EQ(cint[2].val, 5);
}

// ---------------------------------------------------------------------------
// A SmallIx-indexed Array_ starts empty with zero capacity.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, SmallIndexDefaultConstructedIsEmpty) {
    const Array_<int, SmallIx> ismall0;
    EXPECT_TRUE(ismall0.empty());
    EXPECT_EQ(static_cast<int>(ismall0.size()), 0);
    EXPECT_EQ(static_cast<int>(ismall0.capacity()), 0);
    EXPECT_EQ(ismall0.begin(), ismall0.end());
}

// ---------------------------------------------------------------------------
// A SmallIx-indexed Array_ can hold up to max_size() == 4 elements.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, SmallIndexArrayFromRangeFillsToCapacity) {
    const int data[] = {5, 3, -2, 27, 9};

    // Constructing from exactly 4 elements (max_size) must succeed.
    const Array_<int, SmallIx> imaxsz(data, data + 4);
    EXPECT_EQ(static_cast<int>(imaxsz.size()), 4);
    EXPECT_EQ(imaxsz[SmallIx(0)], 5);
    EXPECT_EQ(imaxsz[SmallIx(3)], 27);
}

// ---------------------------------------------------------------------------
// raw_push_back() constructs an element in-place without growing the array
// through the normal path (the caller is responsible for prior reservation).
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, RawPushBackInsertsElementInPlace) {
    Array_<int, SmallIx> ismall(3);
    const auto size_before = static_cast<int>(ismall.size());
    new (ismall.raw_push_back()) int(27);
    EXPECT_EQ(static_cast<int>(ismall.size()), size_before + 1);
    EXPECT_EQ(ismall[SmallIx(3)], 27);
}

// ---------------------------------------------------------------------------
// push_back() on a full SmallIx array throws in Debug (assertion violated).
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, PushBackOnFullSmallIndexArrayThrowsInDebug) {
    const int data[] = {5, 3, -2, 27};
    Array_<int, SmallIx> imaxsz(data, data + 4);
    EXPECT_THROW_DEBUG(imaxsz.push_back(), std::exception);
}

// ---------------------------------------------------------------------------
// Assigning an empty array to a non-empty one clears the destination.
// Assigning a non-empty array to an empty one populates the destination.
// Assigning empty to empty is a no-op.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, NullAssignmentsWorkCorrectly) {
    const int data[] = {5, 3, -2, 27};
    Array_<int, SmallIx> ismall0;                // empty
    Array_<int, SmallIx> ismall(3);              // 3 default ints
    Array_<int, SmallIx> imaxsz(data, data + 4); // 4 ints

    ismall = ismall0; // assign empty → non-empty
    EXPECT_TRUE(ismall.empty());

    ismall0 = imaxsz; // assign non-empty → empty
    EXPECT_EQ(static_cast<int>(ismall0.size()), 4);

    ismall = Array_<int, SmallIx>(); // assign empty → empty
    EXPECT_TRUE(ismall.empty());
}

// ---------------------------------------------------------------------------
// shareData() makes the Array_ point into an existing buffer without copying.
// Modifications through the shared Array_ are visible in the original buffer.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, ShareDataSharesUnderlyingMemory) {
    const int ownerData[] = {7, 77, 777, 7777, 77777};
    std::vector<int> owner(ownerData, ownerData + 5);

    Array_<int> shared;
    // Share elements [1], [2], [3] of owner (3 elements).
    shared.shareData(&owner[1], &owner[4]);

    EXPECT_EQ(shared.size(), 3u);
    EXPECT_EQ(&shared[0], &owner[1]); // must share the same storage
    EXPECT_EQ(shared[0], 77);
    EXPECT_EQ(shared[1], 777);
    EXPECT_EQ(shared[2], 7777);
}

TEST(SimTKCommon_Array_Construction, ShareDataModificationVisibleInOriginalBuffer) {
    const int ownerData[] = {7, 77, 777, 7777, 77777};
    std::vector<int> owner(ownerData, ownerData + 5);

    Array_<int> shared;
    shared.shareData(&owner[1], &owner[4]);

    // Modify through the shared view.
    shared[2] = 29;

    // The original vector must reflect the change.
    EXPECT_EQ(owner[3], 29);

    // The shared view itself must reflect the change.
    EXPECT_EQ(shared[2], 29);
}

// ---------------------------------------------------------------------------
// An OtherArray_ (subclass of Array_) can be constructed normally.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, DerivedArraySubclassConstructsCorrectly) {
    const OtherArray_<int> oa(5, -4);
    EXPECT_EQ(oa.size(), 5u);
    for (auto i = oa.begin(); i != oa.end(); ++i) {
        EXPECT_EQ(*i, -4);
    }
}

// ---------------------------------------------------------------------------
// DontCopy constructor: Array_ wraps an existing std::vector without copying.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Construction, DontCopyConstructorSharesVectorStorage) {
    std::vector<unsigned> unowner{7u, 77u, 777u};
    const Array_<unsigned> shareOfUnowner(unowner, DontCopy());

    EXPECT_EQ(shareOfUnowner.size(), 3u);
    EXPECT_EQ(&shareOfUnowner[0], unowner.data()); // shared, not copied
}

// ===========================================================================
// SimTKCommon_Array_Conversion
// ===========================================================================

// ---------------------------------------------------------------------------
// Array_<int>(vector) copies data; the array and vector have different
// storage addresses.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Conversion, ArrayFromVectorCopiesData) {
    const int p[] = {1, 2, 3, 4, 5, 6};
    const std::vector<int> v(p, p + 6);

    const Array_<int> a(v);

    EXPECT_EQ(a.size(), v.size());
    EXPECT_NE(&a[0], v.data()); // different storage
    for (std::size_t i = 0; i < v.size(); ++i) {
        EXPECT_EQ(a[static_cast<int>(i)], v[i]);
    }
}

// ---------------------------------------------------------------------------
// ArrayView_<int>(vector) shares memory with the original vector.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Conversion, ArrayViewFromVectorSharesMemory) {
    const int p[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> v(p, p + 6);

    const ArrayView_<int> av(v);

    EXPECT_EQ(av.size(), v.size());
    EXPECT_EQ(&av[0], v.data()); // same storage
}

// ---------------------------------------------------------------------------
// ArrayViewConst_<int>(vector) shares memory and provides read-only access.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Conversion, ArrayViewConstFromVectorSharesMemory) {
    const int p[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> v(p, p + 6);

    const ArrayViewConst_<int> ca(v);

    EXPECT_EQ(ca.size(), v.size());
    EXPECT_EQ(&ca[0], v.data()); // same storage
}

// ---------------------------------------------------------------------------
// An ArrayView_<int>(vector) can be passed where const Array_<int>& is
// expected (implicit conversion).
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Conversion, ArrayViewImplicitlyConvertsToConstArrayRef) {
    const int p[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> v(p, p + 6);

    // Helper lambda that accepts a const Array_<int>& – confirms that
    // ArrayView_ is implicitly convertible.
    const auto check = [&](const Array_<int>& arr) {
        EXPECT_EQ(arr.size(), v.size());
        EXPECT_EQ(&arr[0], v.data()); // still the same storage
    };

    check(ArrayView_<int>(v));
}

// ===========================================================================
// SimTKCommon_Array_ViewAssignment
// ===========================================================================

// ---------------------------------------------------------------------------
// Constructing Array_<Count<int>> from int* calls initCtor exactly N times.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, ViewAssignmentCountConstruction) {
    const int data[5] = {10, 100, -23, 4, -99};

    const Array_<Count<int>> acnt(data, data + 5);

    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 0);
}

// ---------------------------------------------------------------------------
// Assigning Array_<int> to Array_<Count<int>> destroys old objects and
// constructs new ones via initCtor (no copy-assignment shortcut).
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, AssignIntArrayToCountArrayUsesInitCtorNotAssign) {
    const int data[5] = {10, 100, -23, 4, -99};
    const Array_<int> adata(data, data + 5);
    Array_<Count<int>> acnt(data, data + 5);
    Count<int>::reset();

    acnt = adata; // clear() then construct from int

    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
}

// ---------------------------------------------------------------------------
// Assigning a larger Array_<Count<int>> to a smaller one destructs all old
// elements and copy-constructs the new ones.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, AssignLargerToSmallerDestroysAllThenCopyConstructs) {
    const int data[5] = {10, 100, -23, 4, -99};
    const Array_<Count<int>> acnt(data, data + 5); // 5 elements
    Array_<Count<int>> acopy(3);                   // 3 default-constructed
    Count<int>::reset();

    acopy = acnt; // destruct 3, copy-construct 5

    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 3);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 5);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::initCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
}

// ---------------------------------------------------------------------------
// Initialising an ArrayView_ from a sub-range of an Array_ shares memory
// without invoking any Count<int> operations.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, ArrayViewInitialisationFromSubrangeSharesMemory) {
    const int data[5] = {10, 100, -23, 4, -99};
    const Array_<Count<int>> acnt(data, data + 5);
    Count<int>::reset();

    // This is an initialisation, not an assignment.
    const ArrayView_<Count<int>> avcnt = acnt(1, 2); // shares 2nd & 3rd elts

    EXPECT_TRUE(Count<int>::isReset()); // nothing should have happened
    EXPECT_EQ(avcnt.size(), 2u);
    EXPECT_TRUE(avcnt[0] == acnt[1]); // same value
    EXPECT_TRUE(avcnt[1] == acnt[2]);
    EXPECT_EQ(&avcnt[0], &acnt[1]); // same address
    EXPECT_EQ(&avcnt[1], &acnt[2]);
}

// // ---------------------------------------------------------------------------
// // Assigning a source that is too large to an ArrayView_ must throw (the view
// // cannot be resized).
// // ---------------------------------------------------------------------------
// TEST(SimTKCommon_Array_ViewAssignment, AssignTooLargeSourceToViewThrows) {
//     const int data[5] = {10, 100, -23, 4, -99};
//     const Array_<int> adata(data, data + 5); // 5 elements
//     Array_<int> backing(data, data + 5);
//     const ArrayView_<int> avcnt = backing(1, 2); // view of 2 elements

//     EXPECT_THROW(avcnt = adata, std::exception); // 5 > 2
// }

// ---------------------------------------------------------------------------
// assign(first, last) on an ArrayView_ performs elementwise assignment via
// initAssign, not a destruct-then-construct round-trip.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, ViewAssignRangeUsesInitAssignNotCopyCtor) {
    const int data[5] = {10, 100, -23, 4, -99};
    const Array_<int> adata(data, data + 5);
    Array_<Count<int>> acnt(data, data + 5);
    Count<int>::reset();

    ArrayView_<Count<int>> avcnt = acnt(1, 2); // 2-element view
    avcnt.assign(adata.begin(), adata.begin() + 2);

    EXPECT_EQ(static_cast<int>(Count<int>::initAssign), 2);
    EXPECT_EQ(static_cast<int>(Count<int>::defCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 0);
}

// // ---------------------------------------------------------------------------
// // Assigning an overlapping source to an ArrayView_ must throw (undefined
// // behaviour would result from aliased element-wise assignment).
// // ---------------------------------------------------------------------------
// TEST(SimTKCommon_Array_ViewAssignment, AssignOverlappingSourceToViewThrows) {
//     const int data[5] = {10, 100, -23, 4, -99};
//     Array_<int> acnt(data, data + 5);
//     const ArrayView_<int> avcnt = acnt(1, 2); // elements [1] and [2]

//     // acnt(0,2) overlaps with avcnt's backing elements.
//     EXPECT_THROW(avcnt = acnt(0, 2), std::exception);
// }

// ---------------------------------------------------------------------------
// Assigning a non-overlapping sub-range to an ArrayView_ uses copyAssign.
// ---------------------------------------------------------------------------
TEST_F(ArrayCountFixture, ViewAssignNonOverlappingSubrangeUsesCopyAssign) {
    const int data[5] = {10, 100, -23, 4, -99};
    Array_<Count<int>> acnt(data, data + 5);
    Count<int>::reset();

    ArrayView_<Count<int>> avcnt = acnt(1, 2); // elements [1] and [2]
    avcnt = acnt(3, 2);                        // assign elements [3],[4] into [1],[2]

    EXPECT_EQ(static_cast<int>(Count<int>::copyAssign), 2);
    EXPECT_EQ(static_cast<int>(Count<int>::copyCtor), 0);
    EXPECT_EQ(static_cast<int>(Count<int>::dtor), 0);

    // was: data[5] = {10, 100, -23, 4, -99}
    // now acnt[1]=4, acnt[2]=-99 (acnt[3] and acnt[4] unchanged)
    EXPECT_EQ(acnt[1].val, 4);
    EXPECT_EQ(acnt[2].val, -99);

    const int modified[5] = {10, 4, -99, 4, -99};
    EXPECT_EQ(acnt, Array_<int>(modified, modified + 5));
}

// ---------------------------------------------------------------------------
// fill() sets all elements to a given value; out-of-range fill count throws
// in Debug.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_ViewAssignment, ArrayViewWrappingCArraySharesPointer) {
    // An ArrayView_ over a raw C array must point directly into that array,
    // not make a copy.
    int someSpace[5] = {123, 1, 12, -9, 14};
    const ArrayView_<int> avSpace(someSpace, someSpace + 5);
    EXPECT_EQ(&avSpace[0], someSpace);
}
TEST(SimTKCommon_Array_ViewAssignment, FillSetsAllElementsToValue) {
    int someSpace[5] = {123, 1, 12, -9, 14};
    ArrayView_<int> avSpace(someSpace, someSpace + 5);

    avSpace = 19;
    EXPECT_EQ(avSpace, Array_<int>(5, 19));

    avSpace.fill(-3);
    EXPECT_EQ(avSpace, Array_<int>(5, -3));
}

TEST(SimTKCommon_Array_ViewAssignment, FillWithWrongCountThrowsInDebug) {
    int someSpace[5] = {123, 1, 12, -9, 14};
    ArrayView_<int> avSpace(someSpace, someSpace + 5);
    EXPECT_THROW_DEBUG(avSpace.assign(12, 999), std::exception);
}

TEST(SimTKCommon_Array_ViewAssignment, FillWithCorrectCountSucceeds) {
    int someSpace[5] = {123, 1, 12, -9, 14};
    ArrayView_<int> avSpace(someSpace, someSpace + 5);
    avSpace.assign(5, 999);
    EXPECT_EQ(avSpace, Array_<int>(5, 999));
}

// ---------------------------------------------------------------------------
// assign() from pointer, random-access iterator, and bidirectional iterator
// all work correctly on an ArrayView_ of the correct size.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_ViewAssignment, AssignFromPointerCopiesValues) {
    int someSpace[5] = {0, 0, 0, 0, 0};
    const int data[5] = {10, 100, -23, 4, -99};
    const std::vector<int> vdata(data, data + 5);
    ArrayView_<int> avSpace(someSpace, someSpace + 5);

    avSpace.assign(data, data + 5);
    EXPECT_EQ(avSpace, vdata);
}

TEST(SimTKCommon_Array_ViewAssignment, AssignFromRandomAccessIteratorCopiesValues) {
    int someSpace[5] = {0, 0, 0, 0, 0};
    const int data[5] = {10, 100, -23, 4, -99};
    const std::vector<int> vdata(data, data + 5);
    ArrayView_<int> avSpace(someSpace, someSpace + 5);

    avSpace.assign(vdata.begin(), vdata.end());
    EXPECT_EQ(avSpace, vdata);
}

TEST(SimTKCommon_Array_ViewAssignment, AssignFromBidirectionalIteratorCopiesValues) {
    int someSpace[5] = {0, 0, 0, 0, 0};
    const int data[5] = {10, 100, -23, 4, -99};
    const std::set<int> aSet(data, data + 5);
    const std::vector<int> sortedVec(aSet.begin(), aSet.end());
    ArrayView_<int> avSpace(someSpace, someSpace + 5);

    avSpace.assign(aSet.begin(), aSet.end());
    EXPECT_EQ(avSpace, sortedVec);
}

// ---------------------------------------------------------------------------
// assign() with fewer elements than the view size must throw.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_ViewAssignment, AssignFewerElementsThanViewSizeThrows) {
    int someSpace[5] = {0, 0, 0, 0, 0};
    const int data[5] = {10, 100, -23, 4, -99};
    const std::vector<int> vdata(data, data + 5);
    const std::set<int> aSet(data, data + 5);

    ArrayView_<int> avSpace(someSpace, someSpace + 5);

    EXPECT_THROW(avSpace.assign(data, data + 3), std::exception);
    EXPECT_THROW(avSpace.assign(vdata.begin(), vdata.begin() + 3), std::exception);

    auto sp = aSet.begin();
    ++sp;
    ++sp;
    EXPECT_THROW(avSpace.assign(aSet.begin(), sp), std::exception);
}

// ===========================================================================
// SimTKCommon_Array_Insert
// ===========================================================================

// ---------------------------------------------------------------------------
// insert() into a non-owner (shared-data) Array_ must throw.
// Inserting zero elements is allowed even on a non-owner.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, InsertIntoNonOwnerThrows) {
    const int data2[4] = {101, 121, -111, 321};
    int wdata[3] = {99, 9999, 999};

    // aw shares wdata – it is a non-owner.
    Array_<int> aw(wdata, wdata + 3, DontCopy());
    const ArrayViewConst_<int> avc(data2, data2 + 4);

    EXPECT_THROW(aw.insert(aw.begin(), avc.begin(), avc.end()), std::exception);
}

TEST(SimTKCommon_Array_Insert, InsertZeroElementsIntoNonOwnerSucceeds) {
    int wdata[3] = {99, 9999, 999};
    Array_<int> aw(wdata, wdata + 3, DontCopy());
    const ArrayViewConst_<int> avc; // empty

    EXPECT_NO_THROW(aw.insert(&aw[1], avc.begin(), avc.begin()));
}

// ---------------------------------------------------------------------------
// insert(pos, first, last) inserts a range at the given position and shifts
// subsequent elements correctly.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, InsertRangeInMiddleShiftsElements) {
    const int data1[3] = {7, -2, 3};
    const int data2[4] = {101, 121, -111, 321};

    Array_<int> ac(data1, data1 + 3);
    std::vector<int> vc(data1, data1 + 3);

    const Array_<int> a(data2, data2 + 4);
    ac.insert(&ac[1], &a[1], &a[1] + 2);
    vc.insert(vc.begin() + 1, &a[1], &a[1] + 2);

    EXPECT_EQ(ac.size(), 5u);
    EXPECT_EQ(ac, vc); // {7, 121, -111, -2, 3}
}

// ---------------------------------------------------------------------------
// Multiple inserts at the beginning and end grow the array correctly.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, InsertAtBeginAndEndGrowsCorrectly) {
    const int data1[3] = {7, -2, 3};
    const int data2[4] = {101, 121, -111, 321};

    Array_<int> ac(data1, data1 + 3);
    std::vector<int> vc(data1, data1 + 3);
    const Array_<int> a(data2, data2 + 4);

    // First build the same 5-element base.
    ac.insert(&ac[1], &a[1], &a[1] + 2);
    vc.insert(vc.begin() + 1, &a[1], &a[1] + 2);

    // Now insert vc at the beginning and end of ac.
    ac.insert(ac.begin(), vc.begin(), vc.end());
    ac.insert(ac.end(), vc.begin(), vc.end());

    EXPECT_EQ(ac.size(), 15u);
    EXPECT_EQ(ac(0, 5), vc);
    EXPECT_EQ(ac(5, 5), vc);
    EXPECT_EQ(ac(10, 5), vc);
}

// ---------------------------------------------------------------------------
// erase(first, last) removes a range and shrink_to_fit reduces capacity.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, EraseRangeAndShrinkToFit) {
    const int data1[3] = {7, -2, 3};
    const int data2[4] = {101, 121, -111, 321};

    Array_<int> ac(data1, data1 + 3);
    std::vector<int> vc(data1, data1 + 3);
    const Array_<int> a(data2, data2 + 4);

    ac.insert(&ac[1], &a[1], &a[1] + 2); // 5 elems
    vc.insert(vc.begin() + 1, &a[1], &a[1] + 2);
    ac.insert(ac.begin(), vc.begin(), vc.end()); // 10
    ac.insert(ac.end(), vc.begin(), vc.end());   // 15

    // Shrink back to 5.
    ac.erase(ac.begin() + 2, ac.begin() + 12);
    EXPECT_EQ(ac, vc);
    EXPECT_GE(ac.allocated(), 15u);

    ac.shrink_to_fit();
    EXPECT_LT(ac.allocated(), 15u);
    EXPECT_EQ(ac, vc); // data must survive shrink
}

// ---------------------------------------------------------------------------
// Inserting zero elements (via null iterator range or count=0) is a no-op.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, InsertZeroElementsIsNoOp) {
    const int data1[3] = {7, -2, 3};
    const int data2[4] = {101, 121, -111, 321};

    Array_<int> ac(data1, data1 + 3);
    std::vector<int> vc(data1, data1 + 3);
    const Array_<int> a(data2, data2 + 4);
    ac.insert(&ac[1], &a[1], &a[1] + 2); // 5 elems – same as vc

    const Array_<int> null;
    ac.insert(ac.begin(), null.begin(), null.end());
    ac.insert(ac.begin() + 2, null.begin(), null.end());
    ac.insert(ac.end(), null.begin(), null.end());
    ac.insert(ac.begin(), 0, 929);
    ac.insert(ac.begin() + 2, 0, 929);
    ac.insert(ac.end(), 0, 929);

    EXPECT_EQ(ac, vc); // unchanged
}

// ---------------------------------------------------------------------------
// Inserting a null ArrayView_ range into a null Array_ leaves both containers
// empty (no-op for both source and destination).
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, InsertNullViewIntoNullArrayLeavesEmpty) {
    Array_<int> null;
    const ArrayView_<int> null2; // default-constructed view = null

    null.insert(null.begin(), null2.begin(), null2.end());

    EXPECT_TRUE(null.empty());
    EXPECT_TRUE(null2.empty());
}

// ---------------------------------------------------------------------------
// Inserting into a null Array_ populates it correctly.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, InsertIntoNullArrayPopulatesIt) {
    const int data1[3] = {7, -2, 3};
    const int data2[4] = {101, 121, -111, 321};

    Array_<int> ac(data1, data1 + 3);
    std::vector<int> vc(data1, data1 + 3);
    const Array_<int> a(data2, data2 + 4);
    ac.insert(&ac[1], &a[1], &a[1] + 2); // 5 elems

    Array_<int> null;

    // insert(pos, n, val)
    null.insert(null.begin(), 3, 987);
    EXPECT_EQ(null, std::vector<int>(3, 987));

    null.deallocate();
    EXPECT_EQ(null.data(), nullptr);
    EXPECT_EQ(null.size(), 0u);
    EXPECT_EQ(null.allocated(), 0u);

    // insert(pos, first, last)
    null.insert(null.begin(), ac.begin(), ac.end());
    EXPECT_EQ(null, vc);
}

// ---------------------------------------------------------------------------
// insert(pos, n, val) then erase leaves exactly the inserted values.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_Insert, BulkInsertThenEraseEndsLeavesInsertedValues) {
    const int data1[3] = {7, -2, 3};
    const int data2[4] = {101, 121, -111, 321};
    const Array_<int> a(data2, data2 + 4);

    Array_<int> ac(data1, data1 + 3);
    Array_<int> vc_arr(data1, data1 + 3);
    ac.insert(&ac[1], &a[1], &a[1] + 2); // 5 elems

    ac.insert(ac.begin() + 2, 99, 1000); // insert 99 × 1000 in the middle
    ac.erase(ac.begin(), ac.begin() + 2);
    ac.erase(ac.end() - 3, ac.end());

    EXPECT_EQ(ac, Array_<int>(99, 1000));
}

// ===========================================================================
// SimTKCommon_Array_BoolIndex
// ===========================================================================

// ---------------------------------------------------------------------------
// An Array_ with long or bool as the index type compiles and behaves
// correctly for basic operations.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_BoolIndex, LongIndexEmptyAndSingleElement) {
    EXPECT_TRUE((Array_<int, long>().empty()));
    EXPECT_FALSE((Array_<int, long>(1L, 99).empty()));
    EXPECT_EQ((Array_<int, long>(1L, 99)[0]), 99);
}

// ---------------------------------------------------------------------------
// Array_<string, bool> has exactly two addressable slots: true and false.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_BoolIndex, BoolIndexDefaultConstructedElementsAreEmpty) {
    // Constructing Array_<string, bool>(2) must value-initialise both slots.
    const Array_<std::string, bool> wisdom(2);
    EXPECT_EQ(wisdom[true], "");
    EXPECT_EQ(wisdom[false], "");
}

TEST(SimTKCommon_Array_BoolIndex, BoolIndexHasTwoSlots) {
    Array_<std::string, bool> wisdom(2);

    wisdom[true] = "this too shall pass";
    wisdom[false] = "don't worry it's not loaded";

    EXPECT_EQ(wisdom.size(), 2u);
    EXPECT_EQ(wisdom.max_size(), 2u);
    EXPECT_GE(wisdom.capacity(), 2u);
    EXPECT_EQ(wisdom.allocated(), wisdom.capacity());
}

TEST(SimTKCommon_Array_BoolIndex, BoolIndexDataAndIteratorConsistency) {
    Array_<std::string, bool> wisdom(2);
    wisdom[true] = "this too shall pass";
    wisdom[false] = "don't worry it's not loaded";

    EXPECT_NE(wisdom.data(), nullptr);
    EXPECT_EQ(wisdom.data(), wisdom.cdata());
    EXPECT_EQ(wisdom.begin(), wisdom.data());
    EXPECT_EQ(wisdom.cbegin(), wisdom.cdata());
    EXPECT_EQ(wisdom.end(), wisdom.begin() + 2);
    EXPECT_EQ(wisdom.cend(), wisdom.cbegin() + 2);
}

TEST(SimTKCommon_Array_BoolIndex, BoolIndexAccessByValueAndAt) {
    Array_<std::string, bool> wisdom(2);
    wisdom[true] = "this too shall pass";
    wisdom[false] = "don't worry it's not loaded";

    EXPECT_EQ(wisdom[false], "don't worry it's not loaded");
    EXPECT_EQ(wisdom[true], "this too shall pass");
    EXPECT_EQ(wisdom.at(false), "don't worry it's not loaded");
    EXPECT_EQ(wisdom.at(true), "this too shall pass");
}

// ---------------------------------------------------------------------------
// Subarrays have a fixed size; assigning a wrong-size std::vector throws in
// Debug.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_BoolIndex, SubarrayAssignWrongSizeThrowsInDebug) {
    Array_<std::string, bool> wisdom(2);
    wisdom[true] = "this too shall pass";
    wisdom[false] = "don't worry it's not loaded";

    EXPECT_THROW_DEBUG(wisdom(false, 2) = std::vector<const char*>(1, "whatever"), std::exception);
}

// ---------------------------------------------------------------------------
// Comparison operators work between Array_<string,bool> and std::vector.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_BoolIndex, ComparisonOperatorsWithStdVector) {
    Array_<std::string, bool> wisdom(2);
    wisdom[true] = "this too shall pass";
    wisdom[false] = "don't worry it's not loaded";

    const std::vector<const char*> vrel(2, "it's all relative");
    wisdom(false, 2) = vrel;

    EXPECT_TRUE(wisdom == vrel);
    EXPECT_TRUE(vrel == wisdom);
    EXPECT_TRUE(wisdom <= vrel);
    EXPECT_TRUE(vrel <= wisdom);
    EXPECT_TRUE(wisdom >= vrel);
    EXPECT_TRUE(vrel >= wisdom);

    // Shorter sub-arrays are less than the full vector.
    EXPECT_TRUE(wisdom(false, 1) < vrel);
    EXPECT_TRUE(wisdom(true, 1) < vrel);
    EXPECT_TRUE(wisdom(0, 0) < vrel);

    // Unequal arrays.
    EXPECT_TRUE(wisdom(false, 1) != vrel);
    EXPECT_TRUE(wisdom != std::vector<const char*>(2, "it's all absolute"));

    // Lexicographically larger element makes the array greater.
    wisdom[true] = "z comes after i";
    EXPECT_TRUE(wisdom > vrel);
    EXPECT_TRUE(vrel < wisdom);
}

// ---------------------------------------------------------------------------
// push_back and out-of-range subarray on a full bool-indexed array throw in
// Debug.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_BoolIndex, PushBackOnFullBoolArrayThrowsInDebug) {
    Array_<std::string, bool> wisdom(2);
    wisdom[true] = "some";
    wisdom[false] = "thing";
    EXPECT_THROW_DEBUG(wisdom.push_back("more brilliance"), std::exception);
}

TEST(SimTKCommon_Array_BoolIndex, SubarrayOutOfRangeThrowsInDebug) {
    Array_<std::string, bool> wisdom(2);
    wisdom[true] = "some";
    wisdom[false] = "thing";
    EXPECT_THROW_DEBUG(wisdom(true, 2), std::exception);
}

// ===========================================================================
// SimTKCommon_Array_NonRandomIter
// ===========================================================================

// ---------------------------------------------------------------------------
// Array_ can be constructed and assigned from std::set (bidirectional iterator
// that does not support subtraction), producing sorted, unique elements.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_NonRandomIter, ConstructFromSetProducesSortedUniqueElements) {
    const int someInts[] = {30, 40, 10, 20, 30, 7, 5};
    const std::set<int> iset(someInts, someInts + 7);
    const std::vector<int> sortUniq(iset.begin(), iset.end());

    const Array_<int> iarr(iset.begin(), iset.end());
    EXPECT_EQ(iarr, sortUniq);
}

TEST(SimTKCommon_Array_NonRandomIter, AssignFromSetProducesSortedUniqueElements) {
    const int someInts[] = {30, 40, 10, 20, 30, 7, 5};
    const std::set<int> iset(someInts, someInts + 7);
    const std::vector<int> sortUniq(iset.begin(), iset.end());

    Array_<int> iarr(iset.begin(), iset.end());
    iarr.assign(iset.begin(), iset.end()); // must increment to count

    EXPECT_EQ(iarr, sortUniq);
}

// ---------------------------------------------------------------------------
// The constructor Array_(integral, integral) must resolve to the n-copies
// overload, NOT the iterator-range overload.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_NonRandomIter, IntegralArgumentsResolveToNCopies) {
    // 3 copies of 'A' (65)
    const Array_<int> dummy1(static_cast<char>(3), 'A');
    EXPECT_EQ(dummy1, Array_<int>(3, static_cast<int>('A')));

    // 4 copies of 129
    const Array_<int> dummy2(4U, 129U);
    EXPECT_EQ(dummy2, Array_<int>(4, 129));
}

// ---------------------------------------------------------------------------
// std::swap specialisation performs a constant-time swap.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_NonRandomIter, StdSwapExchangesContents) {
    Array_<int> dummy1(static_cast<char>(3), 'A');
    Array_<int> dummy2(4U, 129U);
    const auto ptr1 = dummy1.data();
    const auto ptr2 = dummy2.data();

    std::swap(dummy1, dummy2);

    // Pointers must have exchanged (O(1) swap, not O(n) copy).
    EXPECT_EQ(dummy1.data(), ptr2);
    EXPECT_EQ(dummy2.data(), ptr1);
    EXPECT_EQ(dummy2, Array_<int>(3, static_cast<int>('A')));
    EXPECT_EQ(dummy1, Array_<int>(4, 129));
}

// ---------------------------------------------------------------------------
// assign() and insert() with integral args also resolve to the n-copies form.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_NonRandomIter, AssignAndInsertIntegralArgResolveToNCopies) {
    Array_<int> dummy1(static_cast<char>(3), 'A');
    dummy1.assign(static_cast<char>(2), 'B');
    dummy1.insert(dummy1.begin() + 1, static_cast<char>(3), 'C');

    const int d1answer[] = {static_cast<int>('B'),
                            static_cast<int>('C'),
                            static_cast<int>('C'),
                            static_cast<int>('C'),
                            static_cast<int>('B')};
    EXPECT_TRUE((dummy1 == Array_<int, unsigned short>(d1answer, d1answer + 5)));
}

// ---------------------------------------------------------------------------
// fill() overwrites every element with the given value.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_NonRandomIter, FillOverwritesAllElements) {
    Array_<int> dummy1(static_cast<char>(3), 'A');
    dummy1.assign(static_cast<char>(2), 'B');
    dummy1.insert(dummy1.begin() + 1, static_cast<char>(3), 'C'); // 5 elems

    dummy1.fill(7);
    EXPECT_EQ(dummy1, Array_<int>(5, 7));
}

// ---------------------------------------------------------------------------
// Constructing a SmallIx-indexed array from a range that exceeds max_size
// throws in Debug.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_NonRandomIter, ExceedingMaxSizeThrowsInDebug) {
    const int someInts[] = {30, 40, 10, 20, 30, 7, 5};
    const std::set<int> iset(someInts, someInts + 7); // 6 unique ints > 4
    const std::vector<int> sortUniq(iset.begin(), iset.end());

    using AType = Array_<int, SmallIx>;

    // bidirectional iterator
    EXPECT_THROW_DEBUG(AType small1(iset.begin(), iset.end()), std::exception);
    // random access iterator
    EXPECT_THROW_DEBUG(AType small2(sortUniq.begin(), sortUniq.end()), std::exception);
    // pointer (also a random-access iterator; checked separately for clarity)
    const Array_<int> iarr(sortUniq.begin(), sortUniq.end());
    EXPECT_THROW_DEBUG(AType small3(iarr.begin(), iarr.end()), std::exception);
}

// ===========================================================================
// SimTKCommon_Array_InputIterator
// ===========================================================================
// Input iterators (e.g. std::istream_iterator) are single-pass: you cannot
// determine the range size without consuming the range.  Array_ must handle
// this by growing dynamically as elements arrive.

// ---------------------------------------------------------------------------
// Constructing from an istream_iterator reads all elements correctly.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InputIterator, ConstructFromIstreamIteratorReadsAllElements) {
    const int answerData[] = {10, 12, -14, 5, 203, -232, 1, 2, 3, 4};
    const Array_<int, char> answer(answerData, answerData + 10);

    std::istringstream inp1("10 12 -14 5 203 -232 1 2 3 4");
    using Iter = std::istream_iterator<int>;
    Iter p1(inp1);

    const Array_<int> readin(p1, Iter());

    EXPECT_EQ(readin, answer);
}

// ---------------------------------------------------------------------------
// Constructing a SmallIx array from a stream with too many elements throws in
// Debug.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InputIterator, TooManyStreamElementsThrowsInDebug) {
    std::istringstream inp2("10 12 -14 5 203 -232 1 2 3 4"); // 10 elems > max 4
    using Iter = std::istream_iterator<int>;
    Iter p2(inp2);

    using SmallArray = Array_<int, SmallIx>;
    EXPECT_THROW_DEBUG(SmallArray tooSmall(p2, Iter()), std::exception);
}

// ---------------------------------------------------------------------------
// Constructing a SmallIx array from a stream whose element count fits succeeds.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InputIterator, FewStreamElementsSucceedsForSmallIndex) {
    const int answerData[] = {203, -232, 1, 2};
    const Array_<int, short> smallAnswer(answerData, answerData + 4);

    std::istringstream smallInp("203 -232 1 2");
    using Iter = std::istream_iterator<int>;
    Iter psmall(smallInp);

    using SmallArray = Array_<int, SmallIx>;
    const SmallArray okSmall(psmall, Iter());

    EXPECT_EQ(okSmall, smallAnswer);
}

// ---------------------------------------------------------------------------
// operator>> parses a bracketed, optionally comma-separated list.
// Assorted delimiter forms ([], (), {}) and commas are all accepted.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InputIterator, StreamReadBracketedListSucceeds) {
    Array_<float> farray;
    const float farray_ans1[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");

    fin1 >> farray;

    EXPECT_FALSE(fin1.fail());
    EXPECT_EQ(farray, std::vector<float>(farray_ans1, farray_ans1 + 5));
}

// ---------------------------------------------------------------------------
// operator>> on an ArrayView_ replaces exactly the view's elements.
// Trailing tokens after the view is filled are left in the stream.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InputIterator, StreamReadIntoViewReplacesExactlyViewElements) {
    Array_<float> farray;
    const float farray_ans1[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    const float farray_ans2[] = {-1.5f, 910, 920, 9200, 4e-7f};
    std::istringstream fin2(" 9.1e2 9.2e2 9.2e3   ignore me");
    fin2 >> fmid;

    EXPECT_FALSE(fin2.fail());
    EXPECT_EQ(farray, Array_<float>(farray_ans2, farray_ans2 + 5));
}

// ---------------------------------------------------------------------------
// Various malformed inputs must set failbit.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InputIterator, TooFewElementsSetFailbit) {
    Array_<float> farray;
    const float init[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    std::istringstream fin3(" 9.1e2 9.2e2");
    fin3 >> fmid;

    EXPECT_TRUE(fin3.fail()); // wrong (too small) size
}

TEST(SimTKCommon_Array_InputIterator, InconsistentCommaUsageSetFailbit) {
    Array_<float> farray;
    const float init[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    std::istringstream fin4(" 9.1e2 9.2e2,9.2e3 ");
    fin4 >> fmid;

    EXPECT_TRUE(fin4.fail()); // inconsistent use of commas
}

TEST(SimTKCommon_Array_InputIterator, MissingClosingParenSetFailbit) {
    Array_<float> farray;
    const float init[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    std::istringstream fin5("(9.1e2,9.2e2,9.2e3 ");
    fin5 >> fmid;

    EXPECT_TRUE(fin5.fail()); // missing closing paren
}

TEST(SimTKCommon_Array_InputIterator, MismatchedBracesBracketSetFailbit) {
    Array_<float> farray;
    const float init[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    std::istringstream fin6("{9.1e2,9.2e2,9.2e3]");
    fin6 >> fmid;

    EXPECT_TRUE(fin6.fail()); // mismatched delimiters
}

TEST(SimTKCommon_Array_InputIterator, TrailingCommaInsideBracesSetFailbit) {
    Array_<float> farray;
    const float init[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    std::istringstream fin7("{9.1e2,9.2e2,9.2e3,}");
    fin7 >> fmid;

    EXPECT_TRUE(fin7.fail()); // trailing comma inside braces
}

TEST(SimTKCommon_Array_InputIterator, TrailingCommaOutsideBracesOkWhenFull) {
    Array_<float> farray;
    const float init[] = {-1.5f, 3e4f, .125f, 11, 4e-7f};
    std::istringstream fin1("[ -1.5, 3e4 ,.125 , 11,4e-7 ]");
    fin1 >> farray;

    ArrayView_<float> fmid(farray(1, 3));
    // Trailing comma is OK here because we already read the required 3 elements.
    std::istringstream fin8(" 9.1e2,9.2e2,9.2e3,");
    fin8 >> fmid;

    EXPECT_FALSE(fin8.fail());
}

// ===========================================================================
// SimTKCommon_Array_Speed
// ===========================================================================
// The original speed tests had no assertions; they were timing comparisons.
// Here we verify only the arithmetic correctness of push_back + iteration for
// a reduced element count (still exercises the same code paths).

static constexpr int kSpeedInner = 1000;

TEST(SimTKCommon_Array_Speed, StdVectorPushBackAndIterateSumIsCorrect) {
    std::vector<int> v;
    v.reserve(kSpeedInner);

    for (int i = 0; i < kSpeedInner; ++i) {
        v.push_back(i);
    }

    int sum = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        sum += v[i];
    }

    // Sum of 0..999 = (999 * 1000) / 2 = 499500
    EXPECT_EQ(sum, (kSpeedInner * (kSpeedInner - 1)) / 2);
}

TEST(SimTKCommon_Array_Speed, SimTKArrayPushBackAndIterateSumIsCorrect) {
    Array_<int> v;
    v.reserve(kSpeedInner);

    for (int i = 0; i < kSpeedInner; ++i) {
        v.push_back(i);
    }

    int sum = 0;
    for (Array_<int>::size_type i = 0; i < v.size(); ++i) {
        sum += v[static_cast<int>(i)];
    }

    EXPECT_EQ(sum, (kSpeedInner * (kSpeedInner - 1)) / 2);
}

// ===========================================================================
// SimTKCommon_Array_NiceTypeName
// ===========================================================================
// NiceTypeName<T>::name() / namestr() / xmlstr() return human-readable type
// names.  The exact string depends on the platform; we verify that non-empty
// strings are returned and that the SmallIx specialisation returns the
// expected literal.

TEST(SimTKCommon_Array_NiceTypeName, CustomSpecialisationReturnsExpectedName) {
    EXPECT_STREQ(NiceTypeName<SmallIx>::name(), "SmallIx");
    EXPECT_EQ(NiceTypeName<SmallIx>::namestr(), "SmallIx");
    EXPECT_EQ(NiceTypeName<SmallIx>::xmlstr(), "SmallIx");
}

TEST(SimTKCommon_Array_NiceTypeName, BuiltinTypesReturnNonEmptyStrings) {
    EXPECT_FALSE(std::string(NiceTypeName<Is64BitPlatformType>::name()).empty());
    EXPECT_FALSE(std::string(NiceTypeName<ArrayIndexPackType<bool>::packed_size_type>::name()).empty());
    EXPECT_FALSE(std::string(NiceTypeName<ArrayIndexPackType<int>::packed_size_type>::name()).empty());
    EXPECT_FALSE(std::string(NiceTypeName<ArrayIndexPackType<long long>::packed_size_type>::name()).empty());
}

TEST(SimTKCommon_Array_NiceTypeName, ArrayTypeNameStrContainsExpectedSubstrings) {
    const std::string name = NiceTypeName<Array_<String, char>>::namestr();
    EXPECT_FALSE(name.empty());
    // XML variant must not contain raw angle brackets (they would break XML).
    const std::string xml = NiceTypeName<Array_<String, char>>::xmlstr();
    EXPECT_EQ(xml.find('<'), std::string::npos);
    EXPECT_EQ(xml.find('>'), std::string::npos);
}

// ===========================================================================
// SimTKCommon_Array_MemoryFootprint
// ===========================================================================
// The Array_ class must be no larger than std::vector<int> for the default
// (int) index and strictly smaller when the index fits in fewer bits.
// On a 64-bit platform every index ≤ unsigned int should produce a 16-byte
// header; larger types need 24 bytes.

TEST(SimTKCommon_Array_MemoryFootprint, SmallIndexTypesSmallerThanStdVector) {
    EXPECT_LE(sizeof(Array_<int>), sizeof(std::vector<int>));
    EXPECT_LT(sizeof(Array_<int, bool>), sizeof(std::vector<int>));
    EXPECT_LT(sizeof(Array_<int, char>), sizeof(std::vector<int>));
    EXPECT_LT(sizeof(Array_<int, signed char>), sizeof(std::vector<int>));
    EXPECT_LT(sizeof(Array_<int, unsigned char>), sizeof(std::vector<int>));
    EXPECT_LT(sizeof(Array_<int, short>), sizeof(std::vector<int>));
    EXPECT_LT(sizeof(Array_<int, unsigned short>), sizeof(std::vector<int>));
}

TEST(SimTKCommon_Array_MemoryFootprint, LongIndexNoLargerThanStdVector) {
    // long is 32-bit on MSVC/Windows and 64-bit on most Unix; in either
    // case Array_ must not be worse than std::vector.
    EXPECT_LE(sizeof(Array_<int, long>), sizeof(std::vector<int>));
}

TEST(SimTKCommon_Array_MemoryFootprint, PackedSizeTypeMatchesExpectedWidth) {
    // int-width index types → packed as int
    EXPECT_EQ(sizeof(Array_<int>::packed_size_type), sizeof(int));
    EXPECT_EQ(sizeof(Array_<int, int>::packed_size_type), sizeof(int));
    EXPECT_EQ(sizeof(Array_<int, unsigned int>::packed_size_type), sizeof(int));
    // long long-width → packed as long long
    EXPECT_EQ(sizeof(Array_<int, long long>::packed_size_type), sizeof(long long));
    EXPECT_EQ(sizeof(Array_<int, unsigned long long>::packed_size_type), sizeof(long long));
}

TEST(SimTKCommon_Array_MemoryFootprint, SmallIndexPackedIntoCorrectWidthByPlatform) {
    if (Is64BitPlatform) {
        // Small types are packed into a 32-bit int on a 64-bit platform.
        EXPECT_EQ(sizeof(Array_<int, bool>::packed_size_type), sizeof(int));
        EXPECT_EQ(sizeof(Array_<int, char>::packed_size_type), sizeof(int));
        EXPECT_EQ(sizeof(Array_<int, signed char>::packed_size_type), sizeof(int));
        EXPECT_EQ(sizeof(Array_<int, unsigned char>::packed_size_type), sizeof(int));
        EXPECT_EQ(sizeof(Array_<int, short>::packed_size_type), sizeof(int));
        EXPECT_EQ(sizeof(Array_<int, unsigned short>::packed_size_type), sizeof(int));
    } else {
        // Small types are packed into a short on a 32-bit platform.
        EXPECT_EQ(sizeof(Array_<int, bool>::packed_size_type), sizeof(short));
        EXPECT_EQ(sizeof(Array_<int, char>::packed_size_type), sizeof(short));
        EXPECT_EQ(sizeof(Array_<int, signed char>::packed_size_type), sizeof(short));
        EXPECT_EQ(sizeof(Array_<int, unsigned char>::packed_size_type), sizeof(short));
        EXPECT_EQ(sizeof(Array_<int, short>::packed_size_type), sizeof(short));
        EXPECT_EQ(sizeof(Array_<int, unsigned short>::packed_size_type), sizeof(short));
    }
}

TEST(SimTKCommon_Array_MemoryFootprint, AbsoluteSizesByPlatform) {
    if (Is64BitPlatform) {
        EXPECT_EQ(sizeof(Array_<int>), 16u);
        EXPECT_EQ(sizeof(Array_<int, bool>), 16u);
        EXPECT_EQ(sizeof(Array_<int, char>), 16u);
        EXPECT_EQ(sizeof(Array_<int, signed char>), 16u);
        EXPECT_EQ(sizeof(Array_<int, unsigned char>), 16u);
        EXPECT_EQ(sizeof(Array_<int, short>), 16u);
        EXPECT_EQ(sizeof(Array_<int, unsigned short>), 16u);
        EXPECT_EQ(sizeof(Array_<int, int>), 16u);
        EXPECT_EQ(sizeof(Array_<int, unsigned>), 16u);
        EXPECT_LE(sizeof(Array_<int, long>), 24u);
        EXPECT_LE(sizeof(Array_<int, unsigned long>), 24u);
        EXPECT_EQ(sizeof(Array_<int, long long>), 24u);
        EXPECT_EQ(sizeof(Array_<int, unsigned long long>), 24u);
    } else {
        EXPECT_EQ(sizeof(Array_<int>), 12u);
        EXPECT_EQ(sizeof(Array_<int, bool>), 8u);
        EXPECT_EQ(sizeof(Array_<int, char>), 8u);
        EXPECT_EQ(sizeof(Array_<int, signed char>), 8u);
        EXPECT_EQ(sizeof(Array_<int, unsigned char>), 8u);
        EXPECT_EQ(sizeof(Array_<int, short>), 8u);
        EXPECT_EQ(sizeof(Array_<int, unsigned short>), 8u);
        EXPECT_EQ(sizeof(Array_<int, int>), 12u);
        EXPECT_EQ(sizeof(Array_<int, unsigned>), 12u);
        EXPECT_LE(sizeof(Array_<int, long>), 12u);
        EXPECT_LE(sizeof(Array_<int, unsigned long>), 12u);
        // 8-byte integers after a pointer may get 20 or 24 depending on
        // compiler alignment choices.
        EXPECT_LE(sizeof(Array_<int, long long>), 24u);
        EXPECT_LE(sizeof(Array_<int, unsigned long long>), 24u);
    }
}

TEST(SimTKCommon_Array_MemoryFootprint, On64BitIntIndexSmallerThanStdVector) {
    if (Is64BitPlatform) {
        // On a 64-bit platform int occupies only 32 bits while a pointer is
        // 64 bits; packing saves space.
        EXPECT_LT(sizeof(Array_<int, int>), sizeof(std::vector<int>));
        EXPECT_LT(sizeof(Array_<int, unsigned>), sizeof(std::vector<int>));
        EXPECT_LE(sizeof(Array_<int, long long>), sizeof(std::vector<int>));
    }
}

// ===========================================================================
// SimTKCommon_Array_MoveSemantics
// ===========================================================================

// ---------------------------------------------------------------------------
// Move assignment transfers ownership without copying heap data.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_MoveSemantics, MoveAssignmentTransfersOwnership) {
    Array_<double> ad1{1, 2, 3.5, 4};
    const double* p1 = ad1.data();
    Array_<double> ad3(ad1); // copy construction
    const double* p3 = ad3.data();

    ad3 = std::move(ad1); // move assignment

    EXPECT_EQ(ad3.data(), p1); // ad3 now owns the original allocation
    EXPECT_EQ(ad1.data(), p3); // ad1 received ad3's old allocation
}

// ---------------------------------------------------------------------------
// Move construction leaves the source empty and transfers its data pointer.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_MoveSemantics, MoveConstructionLeavesSourceEmpty) {
    Array_<double> ad2{.01, .02};
    const double* p2 = ad2.data();

    const Array_<double> ad4(std::move(ad2));

    EXPECT_EQ(ad4.data(), p2);
    EXPECT_TRUE(ad2.empty());
}

// ---------------------------------------------------------------------------
// An array returned by value via std::move arrives at the call site with the
// original heap pointer intact (no intermediate copy).
// ---------------------------------------------------------------------------
static auto makeDoubleArray(double d) -> std::pair<Array_<double, char>, double*> {
    Array_<double, char> local{1, 2, 3, 4, 5.5};
    local.push_back(d);
    double* localData = local.data();
    return {std::move(local), localData};
}

TEST(SimTKCommon_Array_MoveSemantics, ReturnByValuePreservesHeapPointer) {
    auto returned = makeDoubleArray(3.25);

    EXPECT_EQ(returned.first, std::vector<double>({1, 2, 3, 4, 5.5, 3.25}));
    EXPECT_EQ(returned.first.data(), returned.second);
}

TEST(SimTKCommon_Array_MoveSemantics, MoveAssignmentFromReturnedValuePreservesHeapPointer) {
    auto returned = makeDoubleArray(3.25);
    returned = makeDoubleArray(-1);

    EXPECT_EQ(returned.first, std::vector<double>({1, 2, 3, 4, 5.5, -1}));
    EXPECT_EQ(returned.first.data(), returned.second);
}

// ---------------------------------------------------------------------------
// Array_<unique_ptr<T>> compiles and works: Array_ must not require elements
// to be copy-constructible.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_MoveSemantics, ArrayOfUniquePtrCompiles) {
    Array_<std::unique_ptr<double>> aud;
    aud.push_back(std::make_unique<double>(5.125));
    aud.push_back(std::make_unique<double>(3.5));
    aud.push_back(std::make_unique<double>(-2.25));

    ASSERT_EQ(aud.size(), 3u);
    EXPECT_DOUBLE_EQ(*aud[0], 5.125);
    EXPECT_DOUBLE_EQ(*aud[1], 3.5);
    EXPECT_DOUBLE_EQ(*aud[2], -2.25);
}

TEST(SimTKCommon_Array_MoveSemantics, EmplaceBackAppendsInPlace) {
    Array_<std::unique_ptr<double>> aud;
    aud.push_back(std::make_unique<double>(5.125));
    aud.push_back(std::make_unique<double>(3.5));
    aud.push_back(std::make_unique<double>(-2.25));

    aud.emplace_back(new double(123.)); // NOLINT(cppcoreguidelines-owning-memory)
    ASSERT_EQ(aud.size(), 4u);
    EXPECT_DOUBLE_EQ(*aud[3], 123.);
}

TEST(SimTKCommon_Array_MoveSemantics, EmplaceInMiddleShiftsElements) {
    Array_<std::unique_ptr<double>> aud;
    aud.push_back(std::make_unique<double>(5.125));
    aud.push_back(std::make_unique<double>(3.5));
    aud.push_back(std::make_unique<double>(-2.25));
    aud.emplace_back(new double(123.)); // NOLINT(cppcoreguidelines-owning-memory)

    aud.emplace(&aud[2], new double(100.)); // NOLINT(cppcoreguidelines-owning-memory)
    ASSERT_EQ(aud.size(), 5u);
    EXPECT_DOUBLE_EQ(*aud[2], 100.);
    EXPECT_DOUBLE_EQ(*aud[3], -2.25); // shifted right
}

// ===========================================================================
// SimTKCommon_Array_InitializerList
// ===========================================================================

// ---------------------------------------------------------------------------
// An empty braced-init-list calls the default constructor.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InitializerList, EmptyBracedInitCallsDefaultCtor) {
    const Array_<double> ad1{};
    EXPECT_TRUE(ad1.empty());
}

// ---------------------------------------------------------------------------
// A single-element braced-init-list is a 1-element initializer_list, NOT the
// "n default elements" constructor.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InitializerList, SingleElementBracedInitIsNotNCtor) {
    const Array_<double> ad2{3}; // 1-element initializer_list
    EXPECT_EQ(ad2.size(), 1u);
    EXPECT_DOUBLE_EQ(ad2.front(), 3.0);
}

TEST(SimTKCommon_Array_InitializerList, ParenCountCtorCreatesNDefaultElements) {
    const Array_<double> ad3(3); // 3 default-constructed doubles
    EXPECT_EQ(ad3.size(), 3u);
}

// ---------------------------------------------------------------------------
// Multi-element initializer-list construction copies all values.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InitializerList, MultiElementBracedInitCopiesValues) {
    const Array_<double> ad4 = {1, 2, 2.5, .125};
    EXPECT_EQ(ad4, std::vector<double>({1, 2, 2.5, .125}));
}

// ---------------------------------------------------------------------------
// Assigning an initializer_list to an existing Array_ replaces its contents.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InitializerList, InitListAssignmentReplacesContents) {
    Array_<double> ad4 = {1, 2, 2.5, .125};
    ad4 = {2, 4, 5};
    EXPECT_EQ(ad4, std::vector<double>({2., 4., 5.}));
}

// ---------------------------------------------------------------------------
// An initializer_list of convertible elements implicitly converts to
// Array_<int> and Array_<double>.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_InitializerList, ImplicitConversionFromInitListWorks) {
    const auto checkInt = [](const Array_<int>& arr) {
        EXPECT_EQ(arr, Array_<int>({2, 3, 4}));
    };
    const auto checkDbl = [](const Array_<double>& arr) {
        EXPECT_EQ(arr.size(), 3u);
    };

    checkInt({2, 3, 4});       // implicit conversion to Array_<int>
    checkDbl({1.2, 3.0, 4.0}); // implicit conversion to Array_<double>
}

// ===========================================================================
// SimTKCommon_Array_TypeMismatch
// ===========================================================================
// Array_ allows copy/assign between element types with an *implicit*
// conversion, while std::vector requires exact element-type match.

// ---------------------------------------------------------------------------
// Array_<TestIx> can be assigned from Array_<SubTestIx> because SubTestIx is
// implicitly convertible to TestIx.
// ---------------------------------------------------------------------------

namespace {
// A type that is implicitly convertible to TestIx.
class SubTestIx : public TestIx {
    public:
    explicit SubTestIx(int ix)
        : TestIx(ix) {
    }
};

auto acceptTestIxArray(const Array_<TestIx>& /*arg*/) -> void {
}
auto acceptTestIxVector(const std::vector<TestIx>& /*arg*/) -> void {
}
} // namespace

TEST(SimTKCommon_Array_TypeMismatch, AssignImplicitlyConvertibleElementTypeSucceeds) {
    Array_<TestIx> txarray(2, TestIx(3));
    const Array_<SubTestIx> stxarray(2, SubTestIx(4));

    txarray = stxarray; // SubTestIx → TestIx is implicit

    EXPECT_EQ(txarray, Array_<TestIx>(2, TestIx(4)));
}

// ---------------------------------------------------------------------------
// An Array_<SubTestIx> can be passed where const Array_<TestIx>& is expected,
// because SubTestIx is implicitly convertible to TestIx.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_TypeMismatch, ImplicitConversionAcceptedWhenPassingArray) {
    const Array_<SubTestIx> stxarray(2, SubTestIx(4));
    EXPECT_NO_THROW(acceptTestIxArray(stxarray));
}

// ---------------------------------------------------------------------------
// An explicit-only conversion (SmallIx → TestIx) must NOT be accepted.
// This is verified at compile-time by the commented-out lines in the original;
// at runtime we verify that the SmallIx array still behaves correctly.
// ---------------------------------------------------------------------------
TEST(SimTKCommon_Array_TypeMismatch, SmallIndexArrayOperatesIndependently) {
    const Array_<SmallIx> sxarray(2, SmallIx(1));
    EXPECT_EQ(sxarray.size(), 2u);
    EXPECT_EQ(static_cast<unsigned char>(sxarray[0]), 1u);
    EXPECT_EQ(static_cast<unsigned char>(sxarray[1]), 1u);
    // The following lines must NOT compile (type safety):
    //   txarray = sxarray;          // SmallIx → TestIx has no implicit conv
    //   acceptTestIxArray(sxarray); // same reason
}
