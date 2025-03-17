#!/usr/bin/env python3
"""Classification data structures for the duplicate detection module."""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional
from enum import Enum


class ClassificationType(Enum):
    """Enum representing possible sequence classification types."""
    GOOD = "good"
    BAD = "bad"
    COINFECTION = "coinfection"

# Add decision category enum
class DecisionCategory(Enum):
    """Enum representing possible decision categories in the detection algorithm."""
    SINGLE_SEQUENCE = "SingleSequence"
    BELOW_LOWER_THRESHOLD = "BelowLowerThreshold"
    BELOW_UPPER_THRESHOLD = "BelowUpperThreshold"
    SMALL_CLADE = "SmallClade"
    OUTLIERS_DETECTED = "OutliersDetected"
    TRUE_COINFECTION = "TrueCoinfection"
    NO_TREE = "NoTree"

@dataclass
class Classification:
    """Data class representing a sequence classification result."""
    # Required fields
    sequence_name: str
    classification_type: ClassificationType
    reason: str
    sample_id: str
    group_members: List[str]
    decision_category: DecisionCategory
    contig_stats: Dict[str, float]

    @classmethod
    def from_dict(cls, seq_name: str, data: Dict) -> 'Classification':
        """Create a Classification instance from a dictionary and sequence name.

        Args:
            seq_name: Name of the sequence
            data: Dictionary containing classification information

        Returns:
            Classification object
        """
        # Convert string classification to enum
        classification_str = data.get("classification", "").lower()
        try:
            classification_type = ClassificationType(classification_str)
        except ValueError:
            # Default to BAD classification if invalid
            classification_type = ClassificationType.BAD

        # Extract required fields with defaults
        return cls(
            sequence_name=seq_name,
            classification_type=classification_type,
            reason=data.get("reason", "No reason provided"),
            sample_id=data.get("sample_id", "Unknown"),
            group_members=data.get("group_members", []),
            contig_stats=data.get("contig_stats"),
            decision_category=data.get("decision_category", "")
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert the Classification to a dictionary.

        Returns:
            Dictionary representation of the classification
        """
        result = {
            "classification": self.classification_type.value,
            "reason": self.reason,
            "sample_id": self.sample_id,
            "group_members": self.group_members,
            "decision_category": self.decision_category,
            "contig_stats": self.contig_stats
        }
        return result

    def to_line(self, delimiter: str = "\t") -> str:
        """Format the classification as a delimited line for output files.

        Args:
            delimiter: Field delimiter character (default: tab)

        Returns:
            Delimited string representation of the classification
        """

        # Create line with required fields
        fields = [
            self.sequence_name,
            self.classification_type.value,
            self.decision_category.value,
            str(self.contig_stats) or "NA",
            self.sample_id,
            str(self.group_members),
            self.reason
        ]

        return delimiter.join(fields)

    def to_annotation_line(self, delimiter: str = "\t") -> str:
        """Format the classification as a delimited line for output files.

        Args:
            delimiter: Field delimiter character (default: tab)

        Returns:
            Delimited string representation of the classification
        """

        # Create line with required fields
        fields = [
            self.sequence_name,
            self.classification_type.value,
            self.decision_category.value,
            self.sample_id,
        ]

        return delimiter.join(fields)

    @staticmethod
    def header_line(delimiter: str = "\t") -> str:
        """Get the header line for classification output files.

        Args:
            delimiter: Field delimiter character (default: tab)

        Returns:
            Header line string
        """
        fields = ["sequence_name", "classification", "decision_category",
                "sample_id", "contig_stats", "group_members", "reason",]
        return delimiter.join(fields)

    @staticmethod
    def annotation_line(delimiter: str = "\t") -> str:
        """Get the header line for classification output files.

        Args:
            delimiter: Field delimiter character (default: tab)

        Returns:
            Header line string
        """
        fields = ["sequence_name", "classification", "decision_category", "sample_id"]
        return delimiter.join(fields)

    @property
    def is_good(self) -> bool:
        """Check if this is a 'good' classification."""
        return self.classification_type == ClassificationType.GOOD

    @property
    def is_bad(self) -> bool:
        """Check if this is a 'bad' classification."""
        return self.classification_type == ClassificationType.BAD

    @property
    def is_coinfection(self) -> bool:
        """Check if this is a 'coinfection' classification."""
        return self.classification_type == ClassificationType.COINFECTION

    @property
    def file_classification(self) -> str:
        """Get the file classification ('good' or 'bad') for directory structure.

        Coinfection sequences go in 'good' directory.
        """
        return "good" if self.classification_type in [ClassificationType.GOOD, ClassificationType.COINFECTION] else "bad"
