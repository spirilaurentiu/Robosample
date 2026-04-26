#include <gtest/gtest.h>
#include <string>

#include "SimTKcommon/internal/PrivateImplementation.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include "SimTKcommon.h"

using namespace SimTK;

// =============================================================================
// Forward declarations – mirror what would live in user-visible headers.
// =============================================================================

class Concrete_Impl;
class ConcreteRef_Impl;
class MyHandle_Impl;
class DerivedHandle_Impl;
class DerDerivedHandle_Impl;

// -----------------------------------------------------------------------------
// Concrete – value semantics (third template argument defaults to false)
// -----------------------------------------------------------------------------

/// Handle whose copies are deep (independent) – ordinary object semantics.
class Concrete : public PIMPLHandle<Concrete, Concrete_Impl> {
    public:
    Concrete();

    auto setString(const std::string& n) -> void;
    [[nodiscard]] auto getString() const -> const std::string&;
};

// -----------------------------------------------------------------------------
// ConcreteRef – reference semantics (third template argument == true)
// -----------------------------------------------------------------------------

/// Handle whose copies are shallow (shared) – reference / alias semantics.
class ConcreteRef : public PIMPLHandle<ConcreteRef, ConcreteRef_Impl, true> {
    public:
    ConcreteRef();

    auto setString(const std::string& n) -> void;
    [[nodiscard]] auto getString() const -> const std::string&;
};

// -----------------------------------------------------------------------------
// MyHandle / DerivedHandle / DerDerivedHandle – hierarchical handle classes
// -----------------------------------------------------------------------------

class MyHandle : public PIMPLHandle<MyHandle, MyHandle_Impl, false> {
    public:
    MyHandle();
    explicit MyHandle(MyHandle_Impl* p)
        : HandleBase(p) {
    }

    auto setName(const std::string& n) -> void;
    [[nodiscard]] auto getName() const -> const std::string&;
};

class DerivedHandle : public MyHandle {
    public:
    DerivedHandle();
    explicit DerivedHandle(DerivedHandle_Impl* p);

    auto setReal(Real r) -> void;
    [[nodiscard]] auto getReal() const -> Real;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(DerivedHandle, DerivedHandle_Impl, MyHandle);
};

class DerDerivedHandle : public DerivedHandle {
    public:
    DerDerivedHandle();
    explicit DerDerivedHandle(DerDerivedHandle_Impl* p);

    auto setString(const std::string& s) -> void;
    [[nodiscard]] auto getString() const -> std::string;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(DerDerivedHandle, DerDerivedHandle_Impl, DerivedHandle);
};

// =============================================================================
// Test suite: value-semantic handle (Concrete)
// =============================================================================

// Each handle object should consist of precisely one pointer – no vtable, no
// extra members – regardless of how deep the implementation hierarchy goes.
TEST(SimTKCommon_PrivateImplementation_Concrete, HandleSizeEqualsPointerSize) {
    EXPECT_EQ(sizeof(Concrete), sizeof(void*));
}

// A default-constructed handle must already own a live implementation.
TEST(SimTKCommon_PrivateImplementation_Concrete, DefaultConstructedHandleIsNotEmpty) {
    const Concrete c;
    EXPECT_FALSE(c.isEmptyHandle());
}

// The object that allocated the implementation is its sole owner.
TEST(SimTKCommon_PrivateImplementation_Concrete, DefaultConstructedHandleIsOwner) {
    const Concrete c;
    EXPECT_TRUE(c.isOwnerHandle());
}

// Two independently constructed handles must not alias the same object.
TEST(SimTKCommon_PrivateImplementation_Concrete, TwoDistinctHandlesAreNotTheSameHandle) {
    const Concrete c1;
    const Concrete c2;
    EXPECT_FALSE(c1.isSameHandle(c2));
}

// A handle reports being the same as itself.
TEST(SimTKCommon_PrivateImplementation_Concrete, HandleIsSameAsItself) {
    const Concrete c;
    EXPECT_TRUE(c.isSameHandle(c));
}

// A freshly constructed handle has exactly one handle referencing its impl.
TEST(SimTKCommon_PrivateImplementation_Concrete, InitialImplHandleCountIsOne) {
    const Concrete c;
    EXPECT_EQ(c.getImplHandleCount(), 1);
}

// The stored string is returned correctly after a single set.
TEST(SimTKCommon_PrivateImplementation_Concrete, SetAndGetStringRoundTrips) {
    Concrete c;
    c.setString("hello");
    EXPECT_EQ(c.getString(), "hello");
}

// ------ Copy-construction: value semantics ------

// Copying a Concrete must produce a *distinct* implementation object (deep copy).
TEST(SimTKCommon_PrivateImplementation_Concrete, CopyConstructionProducesDistinctImplementation) {
    Concrete original;
    original.setString("original");

    const Concrete copy(original);

    EXPECT_NE(&original.getImpl(), &copy.getImpl());
}

// Both the original and the copy remain sole owners of their own impls.
TEST(SimTKCommon_PrivateImplementation_Concrete, CopyConstructedHandlesEachHaveImplCountOne) {
    Concrete original;
    const Concrete copy(original);

    EXPECT_EQ(original.getImplHandleCount(), 1);
    EXPECT_EQ(copy.getImplHandleCount(), 1);
}

// The copy starts with the same string value as the original.
TEST(SimTKCommon_PrivateImplementation_Concrete, CopyConstructedHandleHasSameInitialString) {
    Concrete original;
    original.setString("shared_value");

    const Concrete copy(original);

    EXPECT_EQ(copy.getString(), original.getString());
}

// After construction the two handles are fully independent: mutating the copy
// must not change the original.
TEST(SimTKCommon_PrivateImplementation_Concrete, MutatingCopyDoesNotAffectOriginalAfterCopyConstruction) {
    Concrete original;
    original.setString("original");

    Concrete copy(original);
    copy.setString("mutated_copy");

    EXPECT_EQ(original.getString(), "original");
    EXPECT_EQ(copy.getString(), "mutated_copy");
}

// ------ Copy-assignment: value semantics ------

// Copy-assignment must also deep-copy the implementation.
TEST(SimTKCommon_PrivateImplementation_Concrete, CopyAssignmentProducesDistinctImplementation) {
    Concrete source;
    source.setString("source");

    Concrete target;
    target = source;

    EXPECT_NE(&source.getImpl(), &target.getImpl());
}

// Reassigning a handle (clear then assign) still gives an independent impl.
TEST(SimTKCommon_PrivateImplementation_Concrete,
     CopyAssignmentAfterClearHandleProducesDistinctImplementation) {
    Concrete source;
    source.setString("source");

    Concrete target;
    target.clearHandle(); // release previous impl
    target = source;      // deep-copy source's impl

    EXPECT_NE(&source.getImpl(), &target.getImpl());
    EXPECT_EQ(target.getString(), "source");
}

// Mutating the target after assignment must not affect the source.
TEST(SimTKCommon_PrivateImplementation_Concrete, MutatingTargetDoesNotAffectSourceAfterCopyAssignment) {
    Concrete source;
    source.setString("source");

    Concrete target;
    target = source;
    target.setString("target_mutated");

    EXPECT_EQ(source.getString(), "source");
    EXPECT_EQ(target.getString(), "target_mutated");
}

// ------ clearHandle ------

// clearHandle must put the handle into the empty, non-owner state.
TEST(SimTKCommon_PrivateImplementation_Concrete, ClearHandleProducesEmptyNonOwnerHandle) {
    Concrete c;
    c.clearHandle();

    EXPECT_TRUE(c.isEmptyHandle());
    EXPECT_FALSE(c.isOwnerHandle());
}

// ------ Explicit referenceAssign ------

// referenceAssign makes two handles share the same underlying implementation.
// Any mutation through either alias must be visible through the other.
TEST(SimTKCommon_PrivateImplementation_Concrete, ReferenceAssignSharesImplementationAndImplCountIsTwo) {
    Concrete source;
    source.setString("original");

    // Start with an empty handle so we can call referenceAssign.
    Concrete alias;
    alias.clearHandle();
    alias.referenceAssign(source);

    EXPECT_EQ(&alias.getImpl(), &source.getImpl());
    EXPECT_EQ(source.getImplHandleCount(), 2);
}

TEST(SimTKCommon_PrivateImplementation_Concrete, MutationThroughReferenceAssignedAliasIsVisibleOnSource) {
    Concrete source;
    source.setString("before");

    Concrete alias;
    alias.clearHandle();
    alias.referenceAssign(source);
    alias.setString("after");

    EXPECT_EQ(source.getString(), "after");
}

// =============================================================================
// Test suite: reference-semantic handle (ConcreteRef)
// =============================================================================

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, HandleSizeEqualsPointerSize) {
    EXPECT_EQ(sizeof(ConcreteRef), sizeof(void*));
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, DefaultConstructedHandleIsNotEmpty) {
    const ConcreteRef cr;
    EXPECT_FALSE(cr.isEmptyHandle());
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, DefaultConstructedHandleIsOwner) {
    const ConcreteRef cr;
    EXPECT_TRUE(cr.isOwnerHandle());
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, TwoDistinctHandlesAreNotTheSameHandle) {
    const ConcreteRef cr1;
    const ConcreteRef cr2;
    EXPECT_FALSE(cr1.isSameHandle(cr2));
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, HandleIsSameAsItself) {
    const ConcreteRef cr;
    EXPECT_TRUE(cr.isSameHandle(cr));
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, InitialImplHandleCountIsOne) {
    const ConcreteRef cr;
    EXPECT_EQ(cr.getImplHandleCount(), 1);
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, SetAndGetStringRoundTrips) {
    ConcreteRef cr;
    cr.setString("hello");
    EXPECT_EQ(cr.getString(), "hello");
}

// ------ Copy-construction: reference semantics ------

// Copying a ConcreteRef must produce a *shared* implementation (shallow copy).
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, CopyConstructionSharesImplementation) {
    ConcreteRef owner;
    owner.setString("shared");

    const ConcreteRef alias(owner);

    EXPECT_EQ(&owner.getImpl(), &alias.getImpl());
}

// The original is the owner; the copy is a non-owning alias.
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, CopyConstructedAliasIsNotOwner) {
    ConcreteRef owner;
    const ConcreteRef alias(owner);

    EXPECT_TRUE(owner.isOwnerHandle());
    EXPECT_FALSE(alias.isOwnerHandle());
}

// Both handles increment the shared reference count.
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, CopyConstructionIncrementsImplHandleCountToTwo) {
    ConcreteRef owner;
    const ConcreteRef alias(owner);

    EXPECT_EQ(owner.getImplHandleCount(), 2);
    EXPECT_EQ(alias.getImplHandleCount(), 2);
}

// A mutation through one alias is visible through the other – this is the
// defining characteristic of reference semantics.
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, MutationThroughOneAliasIsVisibleThroughAll) {
    ConcreteRef owner;
    owner.setString("original");

    ConcreteRef alias(owner);
    alias.setString("mutated");

    EXPECT_EQ(owner.getString(), "mutated");
}

// ------ Copy-assignment: reference semantics (must clear first) ------

// After clearing, assigning gives a second non-owning alias.
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, CopyAssignmentAfterClearSharesImplementation) {
    ConcreteRef owner;
    owner.setString("shared");

    ConcreteRef alias;
    alias.clearHandle(); // must not be an owner when assigning in ref mode
    alias = owner;

    EXPECT_EQ(&owner.getImpl(), &alias.getImpl());
    EXPECT_FALSE(alias.isOwnerHandle());
}

// Three aliases all point to the same impl and see the same ref-count.
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, ThreeAliasesHaveImplHandleCountOfThree) {
    ConcreteRef owner;
    const ConcreteRef alias1(owner); // copy-construct  → count = 2

    ConcreteRef alias2;
    alias2.clearHandle();
    alias2 = owner; // copy-assign     → count = 3

    EXPECT_EQ(owner.getImplHandleCount(), 3);
    EXPECT_EQ(alias1.getImplHandleCount(), 3);
    EXPECT_EQ(alias2.getImplHandleCount(), 3);
}

// ------ Explicit copyAssign on a reference-semantic handle ------

// copyAssign on a ConcreteRef forces a deep copy (unlike the default
// copy-assign, which makes a new alias).  The resulting handle must own its
// own independent implementation.
TEST(SimTKCommon_PrivateImplementation_ConcreteRef, ExplicitCopyAssignProducesIndependentImplementation) {
    ConcreteRef owner;
    owner.setString("owner_value");

    ConcreteRef deepCopy;
    deepCopy.clearHandle();
    deepCopy.copyAssign(owner); // explicit deep copy

    EXPECT_NE(&owner.getImpl(), &deepCopy.getImpl());
    EXPECT_EQ(deepCopy.getImplHandleCount(), 1);
    EXPECT_EQ(owner.getImplHandleCount(), 1);
}

TEST(SimTKCommon_PrivateImplementation_ConcreteRef, ExplicitCopyAssignDoesNotShareMutations) {
    ConcreteRef owner;
    owner.setString("owner_value");

    ConcreteRef deepCopy;
    deepCopy.clearHandle();
    deepCopy.copyAssign(owner);
    deepCopy.setString("deep_copy_mutated");

    EXPECT_EQ(owner.getString(), "owner_value");
    EXPECT_EQ(deepCopy.getString(), "deep_copy_mutated");
}

// =============================================================================
// Test suite: hierarchical handle classes
// =============================================================================

// All handle types in the hierarchy must still be exactly one pointer wide.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, AllHandleSizesEqualPointerSize) {
    EXPECT_EQ(sizeof(MyHandle), sizeof(void*));
    EXPECT_EQ(sizeof(DerivedHandle), sizeof(void*));
    EXPECT_EQ(sizeof(DerDerivedHandle), sizeof(void*));
}

// A default-constructed base handle is not empty and owns its impl.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DefaultConstructedBaseHandleIsNotEmpty) {
    const MyHandle h;
    EXPECT_FALSE(h.isEmptyHandle());
    EXPECT_TRUE(h.isOwnerHandle());
}

// The base handle stores and retrieves names correctly.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, BaseHandleSetAndGetNameRoundTrips) {
    MyHandle h;
    h.setName("Fred");
    EXPECT_EQ(h.getName(), "Fred");
}

// A derived handle default-constructs without being empty.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DefaultConstructedDerivedHandleIsNotEmpty) {
    const DerivedHandle d;
    EXPECT_FALSE(d.isEmptyHandle());
    EXPECT_TRUE(d.isOwnerHandle());
}

// DerivedHandle stores a Real value (seeded to 27 in its constructor).
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DerivedHandleDefaultRealValueIsTwentySeven) {
    const DerivedHandle d;
    EXPECT_DOUBLE_EQ(d.getReal(), 27.0);
}

// After setting a name on a derived handle the value is retrievable via the
// base-class interface.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DerivedHandleNameIsAccessibleThroughBaseInterface) {
    DerivedHandle d;
    d.setName("Joe");
    EXPECT_EQ(d.getName(), "Joe");
}

// ------ upcast / downcast ------

// Assigning a DerivedHandle to a MyHandle variable (upcast) must produce an
// alias to the same implementation.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, UpcastPreservesImplementationIdentity) {
    DerivedHandle derived;
    derived.setName("Joe");
    derived.setReal(3.14);

    // Copy-construct a base handle from the derived – this is the upcast path.
    const MyHandle base(derived);

    EXPECT_EQ(&base.getImpl(), &derived.getImpl());
}

// Given a MyHandle that was originally a DerivedHandle, downcast must restore
// the derived type without copying.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DowncastOfUpcastedHandleRestoresDerivedType) {
    DerivedHandle original;
    original.setName("Joe");
    original.setReal(99.0);

    const MyHandle asBase(original);                                 // upcast (alias)
    const DerivedHandle& restored = DerivedHandle::downcast(asBase); // downcast (no copy)

    EXPECT_DOUBLE_EQ(restored.getReal(), 99.0);
    EXPECT_EQ(&restored.getImpl(), &original.getImpl());
}

// isInstanceOf must return true when the base handle really holds a derived
// implementation.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, IsInstanceOfReturnsTrueForCorrectDerivedType) {
    DerivedHandle derived;
    const MyHandle asBase(derived);

    EXPECT_TRUE(DerivedHandle::isInstanceOf(asBase));
}

// isInstanceOf must return false when the base handle holds only a base impl.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, IsInstanceOfReturnsFalseForPlainBaseHandle) {
    const MyHandle plain;
    EXPECT_FALSE(DerivedHandle::isInstanceOf(plain));
}

// upcast() on a DerivedHandle must return a reference to the same object
// viewed through the parent type.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, UpcastHelperReturnsSameObjectAsParentReference) {
    DerivedHandle derived;
    derived.setName("upcast_test");

    const MyHandle& viaUpcast = derived.upcast();

    EXPECT_EQ(viaUpcast.getName(), "upcast_test");
    EXPECT_EQ(&viaUpcast.getImpl(), &derived.getImpl());
}

// ------ Two-level derivation ------

// DerDerivedHandle is two levels deep; its own data and all ancestor data must
// be accessible.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DerDerivedHandleDefaultStringIsSet) {
    const DerDerivedHandle dd;
    EXPECT_EQ(dd.getString(), "default dd string");
}

TEST(SimTKCommon_PrivateImplementation_Hierarchy, DerDerivedHandleDefaultRealComesFromDerivedLayer) {
    const DerDerivedHandle dd;
    EXPECT_DOUBLE_EQ(dd.getReal(), 22.345);
}

TEST(SimTKCommon_PrivateImplementation_Hierarchy, DerDerivedHandleSetAndGetString) {
    DerDerivedHandle dd;
    dd.setString("hi there");
    EXPECT_EQ(dd.getString(), "hi there");
}

// Chained upcast: DerDerivedHandle → DerivedHandle → MyHandle.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DoubleUpcastReachesBaseAndPreservesImplementation) {
    DerDerivedHandle dd;
    dd.setName("base_name");

    const MyHandle& atBase = dd.upcast().upcast();

    EXPECT_EQ(atBase.getName(), "base_name");
    EXPECT_EQ(&atBase.getImpl(), &dd.getImpl());
}

// ------ disown / empty handle ------

// disown() transfers ownership to another (empty) handle, leaving the source
// without ownership but still referencing the same impl.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, DisownTransfersOwnershipToTargetHandle) {
    MyHandle owner;
    owner.setName("disown_test");

    MyHandle newOwner(nullptr); // construct an empty handle
    owner.disown(newOwner);

    // The new owner now holds the impl.
    EXPECT_TRUE(newOwner.isOwnerHandle());
    EXPECT_EQ(newOwner.getName(), "disown_test");

    // The original is no longer the owner.
    EXPECT_FALSE(owner.isOwnerHandle());
}

// Constructing a handle with a null pointer produces an empty handle.
TEST(SimTKCommon_PrivateImplementation_Hierarchy, NullPointerConstructorProducesEmptyHandle) {
    const MyHandle empty(nullptr);
    EXPECT_TRUE(empty.isEmptyHandle());
}

// =============================================================================
// Library-side (implementation) code
// =============================================================================
//
// In a real build this would live in a separate translation unit and be linked
// as a binary.  It is included here inline to keep the test self-contained.

// -- Concrete impl -----------------------------------------------------------
class Concrete_Impl : public PIMPLImplementation<Concrete, Concrete_Impl> {
    public:
    [[nodiscard]] auto clone() const -> Concrete_Impl* {
        return new Concrete_Impl(*this);
    }

    std::string s;

    friend class Concrete;
};

namespace SimTK {
template class PIMPLHandle<Concrete, Concrete_Impl>;
template class PIMPLImplementation<Concrete, Concrete_Impl>;
} // namespace SimTK

// -- ConcreteRef impl --------------------------------------------------------
class ConcreteRef_Impl : public PIMPLImplementation<ConcreteRef, ConcreteRef_Impl> {
    public:
    [[nodiscard]] auto clone() const -> ConcreteRef_Impl* {
        return new ConcreteRef_Impl(*this);
    }

    std::string s;

    friend class ConcreteRef;
};

namespace SimTK {
template class PIMPLHandle<ConcreteRef, ConcreteRef_Impl, true>;
template class PIMPLImplementation<ConcreteRef, ConcreteRef_Impl>;
} // namespace SimTK

// -- Concrete & ConcreteRef handle method definitions -----------------------
Concrete::Concrete()
    : HandleBase(new Concrete_Impl()) {
}
ConcreteRef::ConcreteRef()
    : HandleBase(new ConcreteRef_Impl()) {
}

auto Concrete::setString(const std::string& n) -> void {
    updImpl().s = n;
}
auto ConcreteRef::setString(const std::string& n) -> void {
    updImpl().s = n;
}

auto Concrete::getString() const -> const std::string& {
    return getImpl().s;
}
auto ConcreteRef::getString() const -> const std::string& {
    return getImpl().s;
}

// -- MyHandle impl -----------------------------------------------------------
class MyHandle_Impl : public PIMPLImplementation<MyHandle, MyHandle_Impl> {
    public:
    virtual ~MyHandle_Impl() = default;

    [[nodiscard]] virtual auto clone() const -> MyHandle_Impl* {
        return new MyHandle_Impl(*this);
    }

    private:
    std::string name;
    friend class MyHandle;
};

namespace SimTK {
template class PIMPLHandle<MyHandle, MyHandle_Impl>;
template class PIMPLImplementation<MyHandle, MyHandle_Impl>;
} // namespace SimTK

// -- DerivedHandle impl ------------------------------------------------------
class DerivedHandle_Impl : public MyHandle_Impl {
    public:
    [[nodiscard]] auto clone() const -> DerivedHandle_Impl* override {
        return new DerivedHandle_Impl(*this);
    }

    private:
    Real r{};
    friend class DerivedHandle;
};

// -- DerDerivedHandle impl ---------------------------------------------------
class DerDerivedHandle_Impl : public DerivedHandle_Impl {
    public:
    [[nodiscard]] auto clone() const -> DerDerivedHandle_Impl* override {
        return new DerDerivedHandle_Impl(*this);
    }

    private:
    std::string s;
    friend class DerDerivedHandle;
};

// -- MyHandle method definitions ---------------------------------------------
MyHandle::MyHandle()
    : HandleBase(new MyHandle_Impl()) {
}

auto MyHandle::setName(const std::string& n) -> void {
    updImpl().name = n;
}
auto MyHandle::getName() const -> const std::string& {
    return getImpl().name;
}

// -- DerivedHandle method definitions ----------------------------------------
DerivedHandle::DerivedHandle()
    : MyHandle(new DerivedHandle_Impl()) {
    setReal(27);
}

DerivedHandle::DerivedHandle(DerivedHandle_Impl* p)
    : MyHandle(p) {
}

auto DerivedHandle::setReal(Real rr) -> void {
    updImpl().r = rr;
}
auto DerivedHandle::getReal() const -> Real {
    return getImpl().r;
}

// -- DerDerivedHandle method definitions -------------------------------------
DerDerivedHandle::DerDerivedHandle()
    : DerivedHandle(new DerDerivedHandle_Impl()) {
    setString("default dd string");
    setReal(22.345);
}

DerDerivedHandle::DerDerivedHandle(DerDerivedHandle_Impl* p)
    : DerivedHandle(p) {
}

auto DerDerivedHandle::setString(const std::string& s) -> void {
    updImpl().s = s;
}

auto DerDerivedHandle::getString() const -> std::string {
    return getImpl().s;
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(DerivedHandle, DerivedHandle_Impl, MyHandle);
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(DerDerivedHandle, DerDerivedHandle_Impl, DerivedHandle);
