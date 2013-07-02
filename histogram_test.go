package gohistogram

import (
	"math"
	"math/rand"
	"testing"
)

func TestHistogram(t *testing.T) {
	h := NewHistogram(20)
	for i := 0; i < 100; i++ {
		h.Add(rand.NormFloat64())
	}

	if h.total != 100 {
		t.Errorf("Expected h.total to be 100, got ", h.total)
	}

	if per := h.Quantile(0.5); math.Abs(per) > 0.13 {
		t.Errorf("Expected 50th percentile to be 0.0, got %v", per)
	}
	if per := h.Quantile(0.75); math.Abs(per-0.675) > 0.13 {
		t.Errorf("Expected 75th percentile to be 0.675, got %v", per)
	}
	if per := h.Quantile(0.9); math.Abs(per-1.282) > 0.13 {
		t.Errorf("Expected 90th percentile to be 1.282, got %v", per)
	}
}

func TestWeightedHistogram(t *testing.T) {
	h := NewWeightedHistogram(20, 1)
	for i := 0; i < 100; i++ {
		h.Add(rand.NormFloat64())
	}

	if h.total != 100 {
		t.Errorf("Expected h.total to be 100, got ", h.total)
	}

	if per := h.Quantile(0.5); math.Abs(per) > 0.13 {
		t.Errorf("Expected 50th percentile to be 0.0, got %v", per)
	}
	if per := h.Quantile(0.75); math.Abs(per-0.675) > 0.13 {
		t.Errorf("Expected 75th percentile to be 0.675, got %v", per)
	}
	if per := h.Quantile(0.9); math.Abs(per-1.282) > 0.26 {
		t.Errorf("Expected 90th percentile to be 1.282, got %v", per)
	}
}