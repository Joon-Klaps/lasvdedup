import pytest
from lasvdedup.utils.classification import (
    Classification,
    ClassificationType,
    DecisionCategory
)

def test_classification_initialization():
    """Test creating Classification instances with different parameters."""
    # Create a classification with minimal parameters
    classification = Classification(
        sequence_name="seq1",
        classification_type=ClassificationType.GOOD,
        reason="Test reason",
        sample_id="sample1",
        group_members=["seq1", "seq2"],
        decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
        contig_stats={"coverage": 10.5, "identity": 98.2}
    )

    assert classification.sequence_name == "seq1"
    assert classification.classification_type == ClassificationType.GOOD
    assert classification.reason == "Test reason"
    assert classification.sample_id == "sample1"
    assert classification.group_members == ["seq1", "seq2"]
    assert classification.decision_category == DecisionCategory.BELOW_LOWER_THRESHOLD
    assert classification.contig_stats == {"coverage": 10.5, "identity": 98.2}

    # Create a classification with different parameters
    classification = Classification(
        sequence_name="seq2",
        classification_type=ClassificationType.BAD,
        reason="Bad sequence",
        sample_id="sample2",
        group_members=["seq2", "seq3"],
        decision_category=DecisionCategory.BELOW_UPPER_THRESHOLD,
        contig_stats={"coverage": 5.2, "identity": 85.7}
    )

    assert classification.sequence_name == "seq2"
    assert classification.classification_type == ClassificationType.BAD
    assert classification.reason == "Bad sequence"
    assert classification.sample_id == "sample2"
    assert classification.group_members == ["seq2", "seq3"]
    assert classification.decision_category == DecisionCategory.BELOW_UPPER_THRESHOLD
    assert classification.contig_stats == {"coverage": 5.2, "identity": 85.7}

def test_classification_properties():
    """Test Classification helper properties."""
    good = Classification(
        sequence_name="good",
        classification_type=ClassificationType.GOOD,
        reason="Good sequence",
        sample_id="sample1",
        group_members=["good"],
        decision_category=DecisionCategory.SINGLE_SEQUENCE,
        contig_stats={"coverage": 20.0}
    )

    bad = Classification(
        sequence_name="bad",
        classification_type=ClassificationType.BAD,
        reason="Bad sequence",
        sample_id="sample1",
        group_members=["bad"],
        decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
        contig_stats={"coverage": 5.0}
    )

    coinfection = Classification(
        sequence_name="coinf",
        classification_type=ClassificationType.COINFECTION,
        reason="Coinfection",
        sample_id="sample1",
        group_members=["coinf"],
        decision_category=DecisionCategory.TRUE_COINFECTION,
        contig_stats={"coverage": 15.0}
    )

    # Test is_good, is_bad, is_coinfection
    assert good.is_good is True
    assert good.is_bad is False
    assert good.is_coinfection is False

    assert bad.is_good is False
    assert bad.is_bad is True
    assert bad.is_coinfection is False

    assert coinfection.is_good is False
    assert coinfection.is_bad is False
    assert coinfection.is_coinfection is True

    # Test file_classification
    assert good.file_classification == "good"
    assert bad.file_classification == "bad"
    assert coinfection.file_classification == "good"  # Coinfection goes to "good" dir

def test_classification_to_line():
    """Test converting Classification to a delimited line."""
    classification = Classification(
        sequence_name="seq1",
        classification_type=ClassificationType.GOOD,
        reason="Test reason",
        sample_id="sample1",
        group_members=["seq1", "seq2"],
        decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
        contig_stats={"coverage": 10.5}
    )

    line = classification.to_line()

    # Verify line format (tab-delimited)
    parts = line.split("\t")
    assert parts[0] == "seq1"  # sequence_name
    assert parts[1] == "good"  # classification_type.value
    assert parts[2] == "BelowLowerThreshold"  # decision_category.value
    assert parts[3] == "{'coverage': 10.5}" or str(parts[3]) == str({"coverage": 10.5})  # contig_stats
    assert parts[4] == "sample1"  # sample_id
    assert parts[5] == "['seq1', 'seq2']"  # group_members
    assert parts[6] == "Test reason"  # reason

    # Test header line
    header = Classification.header_line()
    assert "sequence_name" in header
    assert "classification" in header
    assert "decision_category" in header
    assert "contig_stats" in header
    assert "sample_id" in header
    assert "group_members" in header
    assert "reason" in header

def test_classification_to_annotation_line():
    """Test converting Classification to an annotation line."""
    classification = Classification(
        sequence_name="seq1",
        classification_type=ClassificationType.GOOD,
        reason="Test reason",
        sample_id="sample1",
        group_members=["seq1", "seq2"],
        decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
        contig_stats={"coverage": 10.5}
    )

    line = classification.to_annotation_line()

    # Verify line format (tab-delimited)
    parts = line.split("\t")
    assert parts[0] == "seq1"  # sequence_name
    assert parts[1] == "good"  # classification_type.value
    assert parts[2] == "BelowLowerThreshold"  # decision_category.value
    assert parts[3] == "sample1"  # sample_id

    # Test annotation header line
    header = Classification.annotation_line()
    assert "sequence_name" in header
    assert "classification" in header
    assert "decision_category" in header
    assert "sample_id" in header

def test_classification_to_dict():
    """Test converting Classification to a dictionary."""
    classification = Classification(
        sequence_name="seq1",
        classification_type=ClassificationType.GOOD,
        reason="Test reason",
        sample_id="sample1",
        group_members=["seq1", "seq2"],
        decision_category=DecisionCategory.BELOW_LOWER_THRESHOLD,
        contig_stats={"coverage": 10.5}
    )

    result = classification.to_dict()

    assert result["classification"] == "good"
    assert result["reason"] == "Test reason"
    assert result["sample_id"] == "sample1"
    assert result["group_members"] == ["seq1", "seq2"]
    assert result["decision_category"] == DecisionCategory.BELOW_LOWER_THRESHOLD
    assert result["contig_stats"] == {"coverage": 10.5}

def test_classification_from_dict():
    """Test creating Classification from a dictionary."""
    data = {
        "classification": "good",
        "reason": "Test reason",
        "sample_id": "sample1",
        "group_members": ["seq1", "seq2"],
        "decision_category": DecisionCategory.BELOW_LOWER_THRESHOLD,
        "contig_stats": {"coverage": 10.5}
    }

    classification = Classification.from_dict("seq1", data)

    assert classification.sequence_name == "seq1"
    assert classification.classification_type == ClassificationType.GOOD
    assert classification.reason == "Test reason"
    assert classification.sample_id == "sample1"
    assert classification.group_members == ["seq1", "seq2"]
    assert classification.decision_category == DecisionCategory.BELOW_LOWER_THRESHOLD
    assert classification.contig_stats == {"coverage": 10.5}

    # Test with invalid classification type
    data["classification"] = "invalid"
    classification = Classification.from_dict("seq1", data)
    assert classification.classification_type == ClassificationType.BAD
