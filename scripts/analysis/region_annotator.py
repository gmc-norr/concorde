"""Region-based variant annotator using BED files and reference FASTA.

Annotates variants with:
- Low complexity region membership (BED lookup)
- GC content in window around variant (reference FASTA)
- Segmental duplication membership (BED lookup)
- Mappability score (BED/bigWig lookup)

All annotations are optional -- if a BED file is not configured, the
corresponding annotation is silently skipped (returns None).
"""

import logging

log = logging.getLogger(__name__)


class RegionAnnotator:
    """Annotate variants with region-based features.

    Args:
        config: Region BED config dict with keys:
            low_complexity, segmental_duplications, mappability
        reference_path: Path to reference FASTA (for GC content)
        gc_window: Window size for GC content calculation (bp)
    """

    def __init__(self, config=None, reference_path=None, gc_window=100):
        self.config = config or {}
        self.reference_path = reference_path
        self.gc_window = gc_window
        self._tabix_cache = {}
        self._fasta = None

    def _get_tabix(self, bed_path):
        """Get or create a pysam TabixFile for the given BED path."""
        if not bed_path:
            return None
        if bed_path not in self._tabix_cache:
            try:
                import pysam

                self._tabix_cache[bed_path] = pysam.TabixFile(str(bed_path))
            except Exception:
                log.warning("Could not open tabix file: %s", bed_path)
                self._tabix_cache[bed_path] = None
        return self._tabix_cache[bed_path]

    def _get_fasta(self):
        """Get or create a pysam FastaFile for GC content."""
        if self._fasta is None and self.reference_path:
            try:
                import pysam

                self._fasta = pysam.FastaFile(str(self.reference_path))
            except Exception:
                log.warning("Could not open reference FASTA: %s", self.reference_path)
                self._fasta = False  # Sentinel to avoid retrying
        return self._fasta if self._fasta is not False else None

    def _tabix_overlaps(self, tabix, chrom, pos):
        """Check if a position overlaps any region in a tabix-indexed BED."""
        if tabix is None:
            return None
        try:
            # BED is 0-based half-open; VCF pos is 1-based
            results = list(tabix.fetch(chrom, pos - 1, pos))
            return len(results) > 0
        except (ValueError, KeyError):
            # Contig not in tabix file
            return None

    def _tabix_score(self, tabix, chrom, pos):
        """Get a numeric score from the 5th column of an overlapping BED region."""
        if tabix is None:
            return None
        try:
            results = list(tabix.fetch(chrom, pos - 1, pos))
            if results:
                fields = results[0].split("\t")
                if len(fields) >= 5:
                    return float(fields[4])
            return None
        except (ValueError, KeyError):
            return None

    def annotate_low_complexity(self, chrom, pos):
        """Check if variant is in a low complexity region."""
        bed_path = self.config.get("low_complexity", "")
        tabix = self._get_tabix(bed_path)
        return self._tabix_overlaps(tabix, chrom, pos)

    def annotate_gc_content(self, chrom, pos, window=None):
        """Compute GC content in a window around the variant position."""
        fasta = self._get_fasta()
        if fasta is None:
            return None

        w = window or self.gc_window
        half = w // 2
        start = max(0, pos - 1 - half)
        end = pos - 1 + half

        try:
            seq = fasta.fetch(chrom, start, end).upper()
            if not seq:
                return None
            gc = sum(1 for b in seq if b in ("G", "C"))
            return gc / len(seq)
        except (ValueError, KeyError):
            return None

    def annotate_segdup(self, chrom, pos):
        """Check if variant is in a segmental duplication region."""
        bed_path = self.config.get("segmental_duplications", "")
        tabix = self._get_tabix(bed_path)
        return self._tabix_overlaps(tabix, chrom, pos)

    def annotate_mappability(self, chrom, pos):
        """Get mappability score for variant position."""
        bed_path = self.config.get("mappability", "")
        tabix = self._get_tabix(bed_path)
        return self._tabix_score(tabix, chrom, pos)

    def annotate_variants(self, variants):
        """Annotate a list of variant objects in-place.

        Sets gc_content, in_low_complexity, in_segdup, mappability_score
        attributes on each variant.

        Args:
            variants: List of variant objects with chrom and pos attributes

        Returns:
            The same list of variants (modified in-place)
        """
        has_lcr = bool(self.config.get("low_complexity", ""))
        has_segdup = bool(self.config.get("segmental_duplications", ""))
        has_mapp = bool(self.config.get("mappability", ""))
        has_ref = bool(self.reference_path)

        if not any([has_lcr, has_segdup, has_mapp, has_ref]):
            log.info("No region annotation sources configured, skipping")
            return variants

        annotated = 0
        for v in variants:
            chrom = getattr(v, "chrom", None)
            pos = getattr(v, "pos", None)
            if chrom is None or pos is None:
                continue

            if has_lcr:
                v.in_low_complexity = self.annotate_low_complexity(chrom, pos)

            if has_ref:
                v.gc_content = self.annotate_gc_content(chrom, pos)

            if has_segdup:
                v.in_segdup = self.annotate_segdup(chrom, pos)

            if has_mapp:
                v.mappability_score = self.annotate_mappability(chrom, pos)

            annotated += 1

        log.info("Annotated %d variants with region features", annotated)
        return variants

    def close(self):
        """Close open file handles."""
        for tbx in self._tabix_cache.values():
            if tbx is not None:
                try:
                    tbx.close()
                except Exception:
                    pass
        if self._fasta and self._fasta is not False:
            try:
                self._fasta.close()
            except Exception:
                pass
