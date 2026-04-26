// TestXml.cpp
//
// Converted from Simbody's SimTK_TEST / SimTK_SUBTEST framework to Google Test.
//
// The original file exercised three distinct subsystems:
//
//   1. convertStringTo<T>()  – generic string-to-type conversion utility
//   2. Xml::Document         – reading and querying XML documents from strings
//   3. Xml::Document         – building XML documents programmatically
//
// Original comments have been preserved wherever possible. Tests that were
// pure I/O (cout only, no assertions) have been augmented with structural and
// value assertions that capture the authors' clear intent.

#include <complex>
#include <cstddef>
#include <gtest/gtest.h>
#include <string>

#include "SimTKcommon/Testing.h"
#include "SimTKcommon/internal/Xml.h"

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

// ============================================================================
// Shared XML document strings
// ============================================================================

// This example is from Wikipedia's XML entry.
// The document intentionally contains duplicate attributes, comments, a CDATA
// section, an unknown tag, and bare top-level text to exercise parser
// robustness.  The top-level prose that precedes <painting> forces the parser
// to wrap the entire content in a synthetic '_Root' element; <painting> is
// therefore a *child* of the root, not the root itself.
static const char* const kXmlPainting = "<?xml version='1.0' encoding='UTF-8'?>\n"
                                        "<!-- a top-level comment -->\n"
                                        "<!-- a multiline\n comment\n third line -->\n"
                                        "     \n"
                                        "but something like this is top level text and will need to get \n"
                                        "moved into a new '_Root' element\n"
                                        "<painting artist='Raphael' artist='metoo'>\n"
                                        "  <img src=\"madonna.jpg\" alt='Foligno Madonna, by Raphael'/>\n"
                                        "  <!-- What follows is a so-called 'caption' -->\n"
                                        "  <caption>This is Raphael's \"Foligno\" Madonna, painted in\n"
                                        "    <date>1511</date>-<date>1512</date>.\n"
                                        "    <![CDATA[some non-Unicode text]]>\n"
                                        "    <  !SOMETHING but this tag is unknown>  \n"
                                        "  </caption>\n"
                                        "  This part is just plain old text.\n"
                                        "  As is this \"quoted\" thing.\n"
                                        "</painting>\n"
                                        "<!whazzis unknown junk>\n"
                                        "<!-- final comment -->";

// A document with no XML elements – only plain prose text.
// The parser should wrap the content in a synthetic root element.
static const char* const kXmlPlainTextFile =
    "That is, the first line should be a declaration, most commonly exactly "
    "the\n"
    "characters shown above, without the \"standalone\" attribute which "
    "will\n"
    "default to \"yes\" anyway. If we don't see a declaration when reading "
    "an XML\n"
    "document, we'll assume we read the one above. Then the document should "
    "contain\n"
    "exactly one top-level (root) element representing the type of document "
    "&amp;\n"
    "document-level attributes.\n";

// A document whose only comment is never closed – must throw on parse.
static const char* const kXmlUnclosedComment =
    "<?xml version='1.0' encoding='UTF-8'?>\n"
    "  <!-- What follows is a so-called 'caption' ->\n" // UNCLOSED!
    "<!whazzis unknown junk>";

// A document consisting of a single comment with no root element.
static const char* const kXmlJustAComment = "<!-- this is the entire contents -->\n";

// White-space-only input – no root element, must throw.
static const char* const kXmlEmpty = "   \n \n \t ";

// ============================================================================
// Internal helpers
// ============================================================================

// Count the number of top-level nodes in a document.
[[nodiscard]] static auto countTopLevelNodes(Xml::Document& doc) -> int {
    int count = 0;
    for (auto it = doc.node_begin(); it != doc.node_end(); ++it) {
        ++count;
    }
    return count;
}

// ============================================================================
// String conversion tests  (original: testStringConvert)
//
// Exercises convertStringTo<T>(), which parses a SimTK::String and returns
// the requested type T.  Covers integers, strings, floats, complex numbers,
// vectors, and fixed/resizable arrays, as well as expected failure modes.
// ============================================================================

// Verify that leading/trailing whitespace is stripped when parsing integers.
TEST(SimTKCommon_StringConvert, IntegerParsedFromStringWithWhitespace) {
    EXPECT_EQ(convertStringTo<int>(" 239\n "), 239);
    EXPECT_EQ(convertStringTo<int>("1234"), 1234);
}

// Verify that String and std::string passthroughs return the original
// content verbatim, including internal whitespace and trailing newlines.
TEST(SimTKCommon_StringConvert, StringPassthroughPreservesContent) {
    EXPECT_EQ(convertStringTo<String>("  lunch box\n"), "  lunch box\n");
    EXPECT_EQ(convertStringTo<std::string>("  lunch box\n"), "  lunch box\n");
}

// Verify that unsigned integers and floats are parsed without data loss.
TEST(SimTKCommon_StringConvert, NumericTypesConvertedCorrectly) {
    EXPECT_EQ(convertStringTo<unsigned>("01234"), 1234U);
    EXPECT_FLOAT_EQ(convertStringTo<float>("1234.5"), 1234.5F);
}

// The parser must reject:
//   • a raw char* target type  (unsupported)
//   • an integer string with trailing non-numeric characters
//   • a floating-point literal aimed at an integer slot
TEST(SimTKCommon_StringConvert, InvalidInputThrows) {
    EXPECT_THROW(convertStringTo<char*>("  lunch box\n"), std::exception);
    EXPECT_THROW(convertStringTo<int>(" 234 j"), std::exception);
    EXPECT_THROW(convertStringTo<int>("345.5"), std::exception);
}

// Verify that a complex number in (re,im) notation is parsed correctly.
TEST(SimTKCommon_StringConvert, ComplexNumberParsed) {
    const std::complex<double> expected{-4.0, 22.0};
    EXPECT_EQ(convertStringTo<std::complex<double>>("(-4,22)"), expected);
}

// Vec3: plain space-separated and comma-separated formats.
TEST(SimTKCommon_StringConvert, Vec3ParsedFromSpaceSeparated) {
    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<Vec3>(\"1 2 3\")",
                                 "Vec3(1,2,3)",
                                 convertStringTo<Vec3>("1 2 3"),
                                 Vec3(1, 2, 3)));

    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<Vec3>(\"1, 2 , 3\")",
                                 "Vec3(1,2,3)",
                                 convertStringTo<Vec3>("1, 2 , 3"),
                                 Vec3(1, 2, 3)));
}

// Vec3: square-bracket forms with and without a leading tilde.
TEST(SimTKCommon_StringConvert, Vec3ParsedFromSquareBrackets) {
    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<Vec3>(\"[ -3 , 5, 6 ] \")",
                                 "Vec3(-3,5,6)",
                                 convertStringTo<Vec3>("[ -3 , 5, 6 ] "),
                                 Vec3(-3, 5, 6)));

    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<Vec3>(\" ~ [ -3 , 5, 6 ] \")",
                                 "Vec3(-3,5,6)",
                                 convertStringTo<Vec3>(" ~ [ -3 , 5, 6 ] "),
                                 Vec3(-3, 5, 6)));
}

// Vec3: round-bracket forms with and without a leading tilde.
TEST(SimTKCommon_StringConvert, Vec3ParsedFromRoundBrackets) {
    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<Vec3>(\"( -3  5 -6 ) \")",
                                 "Vec3(-3,5,-6)",
                                 convertStringTo<Vec3>("( -3  5 -6 ) "),
                                 Vec3(-3, 5, -6)));

    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<Vec3>(\"~( -3  5 -6 ) \")",
                                 "Vec3(-3,5,-6)",
                                 convertStringTo<Vec3>("~( -3  5 -6 ) "),
                                 Vec3(-3, 5, -6)));
}

// Mismatched opening/closing brackets must throw.
TEST(SimTKCommon_StringConvert, Vec3ThrowsOnMismatchedBrackets) {
    EXPECT_THROW(convertStringTo<Vec3>("( -3  5 -6 ] "), std::exception);
    EXPECT_THROW(convertStringTo<Vec3>(" -3  5 -6 ] "), std::exception);
    // Bare '~' with no bracket is also illegal.
    EXPECT_THROW(convertStringTo<Vec3>(" ~ -3  5 -6 "), std::exception);
}

// Verify that a Vec<2, complex<float>> can be parsed from bracket notation.
TEST(SimTKCommon_StringConvert, ComplexVec2Parsed) {
    using fCVec2 = Vec<2, std::complex<float>>;
    const fCVec2 expected{std::complex<float>(1.F, 2.F), std::complex<float>(3.F, 4.F)};
    EXPECT_TRUE(AssertSimTKEqual("convertStringTo<fCVec2>(\"[(1,2) (3,4)]\")",
                                 "fCVec2((1,2),(3,4))",
                                 convertStringTo<fCVec2>("[(1,2) (3,4)]"),
                                 expected));
}

// Verify that a space-separated list of integers is parsed into an Array_.
TEST(SimTKCommon_StringConvert, ArrayIntParsed) {
    const Array_<int> result = convertStringTo<Array_<int>>("1 2 3 4");
    ASSERT_EQ(result.size(), 4U);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 2);
    EXPECT_EQ(result[2], 3);
    EXPECT_EQ(result[3], 4);
}

// An ArrayView_ has fixed capacity; assigning more elements than it can hold
// must throw because ArrayView_ is fixed size.
TEST(SimTKCommon_StringConvert, ArrayViewThrowsWhenSourceIsTooLarge) {
    Array_<float> af(2);
    EXPECT_THROW(String(" -.25, .5, 29.2e4 ").convertTo<ArrayView_<float>>(af), std::exception);
}

// An Array_ is resizable; assigning more elements than the current size must
// succeed and the container must hold the correct values afterwards.
TEST(SimTKCommon_StringConvert, ResizableArrayAcceptsLargerSource) {
    Array_<float> af(2);
    ASSERT_NO_THROW(String(" -.25, .5, 29.2e4 ").convertTo<Array_<float>>(af));
    ASSERT_EQ(af.size(), 3U);
    EXPECT_FLOAT_EQ(af[0], -.25F);
    EXPECT_FLOAT_EQ(af[1], .5F);
    EXPECT_FLOAT_EQ(af[2], 292000.0F);
}

// ============================================================================
// XML document tests – reading from strings  (original: testXmlFromString)
//
// Exercises Xml::Document parsing from string inputs of varying validity:
// comment-only, plain-text, well-formed, empty, and malformed documents.
// Also covers node/element iteration, attribute access, and serialisation.
// ============================================================================

// A document that contains only a comment (no root element) must be readable
// without throwing.
TEST(SimTKCommon_XmlFromString, ParsesCommentOnlyDocument) {
    Xml::Document doc;
    EXPECT_NO_THROW(doc.readFromString(kXmlJustAComment));
}

// Parsing a plain-text file (no XML tags at all) must succeed; the parser
// wraps the content in a synthetic root element.
TEST(SimTKCommon_XmlFromString, ParsesPlainTextFile) {
    Xml::Document doc;
    EXPECT_NO_THROW(doc.readFromString(kXmlPlainTextFile));
}

// Note that the "condense white space" setting is global, not
// document-specific.  Toggling it must be immediately reflected by the query
// function.
TEST(SimTKCommon_XmlFromString, WhiteSpaceCondensationIsGlobalSetting) {
    Xml::Document::setXmlCondenseWhiteSpace(false);
    EXPECT_FALSE(Xml::Document::isXmlWhiteSpaceCondensed());

    // Restore the default so subsequent tests are unaffected.
    Xml::Document::setXmlCondenseWhiteSpace(true);
    EXPECT_TRUE(Xml::Document::isXmlWhiteSpaceCondensed());
}

// White-space-only input has no root element and must throw.
TEST(SimTKCommon_XmlFromString, ThrowsOnEmptyDocument) {
    Xml::Document doc;
    EXPECT_THROW(doc.readFromString(kXmlEmpty), std::exception);
}

// A document with an unclosed comment must throw on parse.
TEST(SimTKCommon_XmlFromString, ThrowsOnUnclosedComment) {
    Xml::Document doc;
    EXPECT_THROW(doc.readFromString(kXmlUnclosedComment), std::exception);
}

// After parsing the painting document, the XML prolog values (version and
// encoding) must match those declared in the header.
TEST(SimTKCommon_XmlFromString, PaintingDocumentHasCorrectXmlProlog) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    EXPECT_EQ(doc.getXmlVersion(), "1.0");
    EXPECT_EQ(doc.getXmlEncoding(), "UTF-8");
    // getXmlIsStandalone() must be callable without throwing even though no
    // explicit standalone attribute was declared in the prolog.
    EXPECT_NO_THROW(doc.getXmlIsStandalone());
}

// The root element of the painting document must have at least one child node
// and at least one comment child.  Because top-level prose precedes <painting>
// the parser wraps everything in a synthetic '_Root'; hasNode() applies to
// whatever getRootElement() returns.
TEST(SimTKCommon_XmlFromString, PaintingRootElementHasChildNodes) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    const Xml::Element root = doc.getRootElement();
    EXPECT_TRUE(root.hasNode());
    EXPECT_TRUE(root.hasNode(Xml::CommentNode));
}

// The <painting> element (child of the synthetic '_Root') must have exactly
// two child elements: <img> and <caption>, in that order.
TEST(SimTKCommon_XmlFromString, PaintingElementHasImgAndCaptionChildren) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    // Top-level text before <painting> causes a '_Root' wrapper; <painting>
    // is therefore a child of the root, not the root itself.
    Xml::Element root = doc.getRootElement();
    Xml::Element painting = root.getRequiredElement("painting");
    const Array_<Xml::Element> children = painting.getAllElements();

    ASSERT_EQ(children.size(), 2U);
    EXPECT_EQ(children[0].getElementTag(), "img");
    EXPECT_EQ(children[1].getElementTag(), "caption");
}

// element_begin("caption") on the painting element must yield exactly one hit.
TEST(SimTKCommon_XmlFromString, PaintingHasExactlyOneCaptionElement) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    Xml::Element painting = doc.getRootElement().getRequiredElement("painting");

    int count = 0;
    for (auto ep = painting.element_begin("caption"); ep != painting.element_end(); ++ep) {
        ++count;
    }
    EXPECT_EQ(count, 1);
}

// The <img> element inside <painting> must expose a src attribute with the
// correct value, and an alt attribute.
TEST(SimTKCommon_XmlFromString, PaintingImgHasCorrectAttributes) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    Xml::Element img = doc.getRootElement().getRequiredElement("painting").getRequiredElement("img");

    const Array_<Xml::Attribute> attrs = img.getAllAttributes();
    ASSERT_GE(attrs.size(), 2U);

    bool hasSrc = false;
    bool hasAlt = false;
    for (const auto& attr : attrs) {
        if (attr.getName() == "src") {
            hasSrc = true;
            EXPECT_EQ(attr.getValue(), "madonna.jpg");
        }
        if (attr.getName() == "alt") {
            hasAlt = true;
        }
    }
    EXPECT_TRUE(hasSrc);
    EXPECT_TRUE(hasAlt);
}

// A compact serialisation must be strictly shorter than a pretty-printed one.
TEST(SimTKCommon_XmlFromString, CompactStringIsShorterThanPrettyString) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    String pretty;
    String compact;
    doc.writeToString(pretty);
    doc.writeToString(compact, /*compact=*/true);

    EXPECT_LT(compact.size(), pretty.size());
}

// getRequiredElement() must throw when the requested child tag does not exist.
TEST(SimTKCommon_XmlFromString, GetRequiredElementThrowsWhenMissing) {
    Xml::Document doc;
    ASSERT_NO_THROW(doc.readFromString(kXmlPainting));

    Xml::Element root = doc.getRootElement();
    EXPECT_THROW(root.getRequiredElement("nonexistent"), std::exception);
}

// ============================================================================
// XML document tests – construction from scratch  (original: testXmlFromScratch)
//
// Exercises programmatic creation of Xml::Documents and elements, including
// orphan management, element value and attribute manipulation, node
// insertion/removal ordering, deep copy independence, and element cloning.
// ============================================================================

// A freshly created document must report exactly the root tag that was set.
TEST(SimTKCommon_XML, CreatesDocumentWithCustomRootTag) {
    Xml::Document doc;
    doc.setRootTag("MyDoc");
    EXPECT_EQ(doc.getRootTag(), "MyDoc");
}

// Nodes that have not been inserted into any document are orphans; once
// inserted via insertTopLevelNodeAfter they must no longer be orphans.
TEST(SimTKCommon_XML, OrphanStateChangesAfterInsertion) {
    Xml::Document doc;
    doc.setRootTag("Root");

    Xml::Comment c("This is a comment.");
    EXPECT_TRUE(c.isOrphan());

    doc.insertTopLevelNodeAfter(doc.node_begin(), c);
    EXPECT_FALSE(c.isOrphan());
}

// A newly created element has an empty value; updValue() appends in place
// while setValue() replaces any previous content entirely.
TEST(SimTKCommon_XML, ElementValueCanBeSetAndAppended) {
    Xml::Element e("elementTag");

    // A newly created element has an empty value.
    EXPECT_EQ(e.getValue(), "");

    e.updValue() += "AVALUE:";
    EXPECT_EQ(e.getValue(), "AVALUE:");

    e.setValue("this is the only value");
    e.updValue() += " (but then I added this)";
    EXPECT_EQ(e.getValue(), "this is the only value (but then I added this)");

    // setValue() replaces any previous content.
    e.setValue("9 10 -3.2e-4");
    EXPECT_EQ(e.getValue(), "9 10 -3.2e-4");

    // We're never going to insert this element so its heap space will
    // leak if we don't explicitly call clearOrphan().
    e.clearOrphan();
}

// An attribute stored via setAttributeValue as a String representation of a
// Vec2 must be retrievable via getRequiredAttributeValueAs<Vec2>() and
// compare equal to the original value.
TEST(SimTKCommon_XML, AttributeRoundTripsAsVec2) {
    Xml::Element e("elementTag");

    e.setAttributeValue("attr1", String(Vec2(9, -9)));
    const Vec2 retrieved = e.getRequiredAttributeValueAs<Vec2>("attr1");

    EXPECT_TRUE(AssertSimTKEqual("getRequiredAttributeValueAs<Vec2>(\"attr1\")",
                                 "Vec2(9,-9)",
                                 retrieved,
                                 Vec2(9, -9)));

    e.clearOrphan();
}

// Inserting a top-level node must increase the document's top-level node
// count by exactly one.
TEST(SimTKCommon_XML, NodeInsertionIncreasesTopLevelNodeCount) {
    Xml::Document doc;
    doc.setRootTag("Root");

    const int before = countTopLevelNodes(doc);

    Xml::Comment c("This is a comment.");
    doc.insertTopLevelNodeAfter(doc.node_begin(), c);

    EXPECT_EQ(countTopLevelNodes(doc), (before + 1));
}

// Deep copy via operator= must be independent of the original: erasing a node
// from the original must leave the copy's node count unchanged.
TEST(SimTKCommon_XML, DeepCopyIsIndependentOfOriginal) {
    Xml::Document original;
    original.setRootTag("Root");

    Xml::Comment c("A comment.");
    original.insertTopLevelNodeAfter(original.node_begin(), c);

    const int originalCount = countTopLevelNodes(original);

    Xml::Document copy = original; // deep copy via assignment operator

    original.eraseTopLevelNode(original.node_begin());

    EXPECT_EQ(countTopLevelNodes(original), (originalCount - 1))
        << "Original should have one fewer node after erasure.";
    EXPECT_EQ(countTopLevelNodes(copy), originalCount)
        << "Deep copy must retain all nodes independently of the original.";
}

// eraseTopLevelNode() must reduce the document's top-level node count by
// exactly one.
TEST(SimTKCommon_XML, EraseTopLevelNodeReducesNodeCount) {
    Xml::Document doc;
    doc.setRootTag("Root");

    Xml::Comment c("Comment to erase.");
    doc.insertTopLevelNodeAfter(doc.node_begin(), c);

    const int before = countTopLevelNodes(doc);
    doc.eraseTopLevelNode(doc.node_begin());

    EXPECT_EQ(countTopLevelNodes(doc), (before - 1));
}

// An element inserted into another element via insertNodeAfter must be
// findable afterwards via getRequiredElement().
TEST(SimTKCommon_XML, ChildElementIsFoundAfterInsertion) {
    Xml::Document doc;
    doc.setRootTag("Root");
    Xml::Element root = doc.getRootElement();

    Xml::Element child("childTag");
    root.insertNodeAfter(root.node_begin(), child);

    EXPECT_NO_THROW(root.getRequiredElement("childTag"));
}

// An element initialised with a Vec3 value and then inserted as a child must
// round-trip that value correctly through getValueAs<Vec3>().
TEST(SimTKCommon_XML, Vec3ElementValuePreservedAfterInsertion) {
    Xml::Document doc;
    doc.setRootTag("Root");
    Xml::Element root = doc.getRootElement();

    Xml::Element parent("parent");
    root.insertNodeAfter(root.node_begin(), parent);

    Xml::Element child("anotherElt", Vec3(.1, .2, .3));
    parent.insertNodeAfter(parent.element_end(), child);

    const Xml::Element found = parent.getRequiredElement("anotherElt");
    EXPECT_TRUE(AssertSimTKEqual("anotherElt getValueAs<Vec3>()",
                                 "Vec3(0.1, 0.2, 0.3)",
                                 found.getValueAs<Vec3>(),
                                 Vec3(.1, .2, .3)));
}

// clone() must produce an element with the same tag that is structurally
// independent of the original: mutating the clone's tag must not affect the
// original, and vice versa.
TEST(SimTKCommon_XML, CloneIsIndependentOfOriginal) {
    Xml::Element original("originalTag");
    Xml::Element cloned = original.clone();
    cloned.setElementTag("clonedTag");

    EXPECT_EQ(original.getElementTag(), "originalTag");
    EXPECT_EQ(cloned.getElementTag(), "clonedTag");

    original.clearOrphan();
    cloned.clearOrphan();
}

// removeNode() must detach the named element from the tree so that a
// subsequent getRequiredElement() call for the same tag throws.
TEST(SimTKCommon_XML, RemoveNodeExtractsElement) {
    Xml::Document doc;
    doc.setRootTag("Root");
    Xml::Element root = doc.getRootElement();

    Xml::Element e("targetElement");
    root.insertNodeAfter(root.node_begin(), e);

    // Verify the element is present before removal.
    ASSERT_NO_THROW(root.getRequiredElement("targetElement"));

    Xml::Node extracted = root.removeNode(root.element_begin("targetElement"));

    // After removal the element must no longer be found.
    EXPECT_THROW(root.getRequiredElement("targetElement"), std::exception);

    // Free the extracted orphan to avoid a memory leak (see original comment
    // about neverMind.clearOrphan()).
    extracted.clearOrphan();
}