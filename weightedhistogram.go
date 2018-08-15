// Package gohistogram contains implementations of weighted and exponential histograms.
package gohistogram

// Copyright (c) 2013 VividCortex, Inc. All rights reserved.
// Please see the LICENSE file for applicable license terms.

import (
	"encoding/json"
	"fmt"

	"sort"

	"github.com/ventu-io/slf"
)

var log = slf.WithContext("gohistogram")

// A WeightedHistogram implements Histogram. A WeightedHistogram has bins that have values
// which are exponentially weighted moving averages. This allows you keep inserting large
// amounts of data into the histogram and approximate quantiles with recency factored in.
type WeightedHistogram struct {
	bins    []bin
	maxbins int
	total   int64
	alpha   float64
}

type histogramStateBin struct {
	Value float64
	Count float64
}

type histogramState struct {
	Bins    []histogramStateBin
	MaxBins int
	Total   int64
	Alpha   float64
}

//MarshalJSON implements json.Marshaller
func (h *WeightedHistogram) MarshalJSON() ([]byte, error) {
	hs := &histogramState{
		Bins:    make([]histogramStateBin, len(h.bins), len(h.bins)),
		MaxBins: h.maxbins,
		Total:   h.total,
		Alpha:   h.alpha,
	}
	for i, v := range h.bins {
		hs.Bins[i] = histogramStateBin{
			Value: v.value,
			Count: v.count,
		}
	}
	result, err := json.Marshal(hs)
	if err != nil {
		log.Errorf("MarshalJSON: %+v error %s", hs, err)
		return nil, err
	}
	return result, nil
}

//UnmarshalJSON implements json.Unmarshaller
func (h *WeightedHistogram) UnmarshalJSON(data []byte) error {
	var hs histogramState
	err := json.Unmarshal(data, &hs)
	if err != nil {
		return err
	}

	h.bins = make([]bin, len(hs.Bins), len(hs.Bins))
	h.maxbins = hs.MaxBins
	h.alpha = hs.Alpha
	h.total = hs.Total

	for i, v := range hs.Bins {
		h.bins[i] = bin{
			value: v.Value,
			count: v.Count,
		}
	}

	return nil
}

// NewWeightedHistogram returns a new WeightedHistogram with a maximum of n bins with a decay factor
// of alpha.
//
// There is no "optimal" bin count, but somewhere between 20 and 80 bins should be
// sufficient.
//
// Alpha should be set to 1 - (2 / (N+1)), where N represents the average age of the moving window.
// For example, a 60-second window with an average age of 30 seconds would yield an
// alpha of 0.935483870967742.
func NewWeightedHistogram(n int, alpha float64) *WeightedHistogram {
	return &WeightedHistogram{
		bins:    make([]bin, 0),
		maxbins: n,
		total:   0,
		alpha:   alpha,
	}
}

func ewma(existingVal float64, newVal float64, alpha float64) (result float64) {
	result = newVal*(1-alpha) + existingVal*alpha
	return
}

func (h *WeightedHistogram) scaleDown(except int) {
	for i := range h.bins {
		if i != except {
			h.bins[i].count = ewma(h.bins[i].count, 0, h.alpha)
		}
	}
}

//Add adds value to histogram
func (h *WeightedHistogram) Add(n float64) {
	defer h.trim()
	for i := range h.bins {
		if h.bins[i].value == n {
			h.bins[i].count++

			defer h.scaleDown(i)
			return
		}

		if h.bins[i].value > n {

			newbin := bin{value: n, count: 1}
			head := append(make([]bin, 0), h.bins[0:i]...)

			head = append(head, newbin)
			tail := h.bins[i:]
			h.bins = append(head, tail...)

			defer h.scaleDown(i)
			return
		}
	}

	h.bins = append(h.bins, bin{count: 1, value: n})
}

// Quantile implements Histogram.Quantile and returns an approximation.
func (h *WeightedHistogram) Quantile(q float64) float64 {
	count := q * float64(h.total)
	for i := range h.bins {
		count -= float64(h.bins[i].count)

		if count <= 0 {
			return h.bins[i].value
		}
	}

	return -1
}

// CDF returns the value of the cumulative distribution function
// at x
func (h *WeightedHistogram) CDF(x float64) float64 {
	count := 0.0
	for i := range h.bins {
		if h.bins[i].value <= x {
			count += float64(h.bins[i].count)
		}
	}

	return count / float64(h.total)
}

// Mean returns the sample mean of the distribution
func (h *WeightedHistogram) Mean() float64 {
	if h.total == 0 {
		return 0
	}

	sum := 0.0

	for i := range h.bins {
		sum += h.bins[i].value * h.bins[i].count
	}

	return sum / float64(h.total)
}

// Modes returns values for first n maximums from histogram
func (h *WeightedHistogram) Modes(n int) []float64 {
	result := make([]float64, 0)
	if h.total == 0 {
		return result
	}
	tmp := make([]bin, 0)
	tmp = append(tmp, h.bins...)
	sort.Slice(tmp, func(i, j int) bool {
		return tmp[i].count >= tmp[j].count
	})
	for i := 0; i < n && i < len(tmp); i++ {
		result = append(result, tmp[i].value)
	}
	return result
}

// Variance returns the variance of the distribution
func (h *WeightedHistogram) Variance() float64 {
	if h.total == 0 {
		return 0
	}

	sum := 0.0
	mean := h.Mean()

	for i := range h.bins {
		sum += (h.bins[i].count * (h.bins[i].value - mean) * (h.bins[i].value - mean))
	}

	return sum / float64(h.total)
}

// Count returns approximate decayed data count
func (h *WeightedHistogram) Count() int64 {
	return h.total
}

func (h *WeightedHistogram) trim() {
	total := 0.0
	for i := range h.bins {
		total += h.bins[i].count
	}
	h.total = int64(total)
	for len(h.bins) > h.maxbins {

		// Find closest bins in terms of value
		minDelta := 1e99
		minDeltaIndex := 0
		for i := range h.bins {
			if i == 0 {
				continue
			}

			if delta := h.bins[i].value - h.bins[i-1].value; delta < minDelta {
				minDelta = delta
				minDeltaIndex = i
			}
		}

		// We need to merge bins minDeltaIndex-1 and minDeltaIndex
		b1 := h.bins[minDeltaIndex-1]
		b2 := h.bins[minDeltaIndex]
		totalCount := b1.count + b2.count
		var newValue float64
		if totalCount <= 1 {
			newValue = (b1.value + b2.value) / 2
		} else {
			newValue = (b1.value*b1.count + b2.value*b2.count) / totalCount // weighted average
		}
		//log.Debugf("trim: %d %d %f %f", len(h.bins), minDeltaIndex, newValue, totalCount)
		mergedbin := bin{
			value: newValue,
			count: totalCount, // summed heights
		}
		head := append(make([]bin, 0), h.bins[0:minDeltaIndex-1]...)
		tail := append([]bin{mergedbin}, h.bins[minDeltaIndex+1:]...)
		h.bins = append(head, tail...)
	}
}

// String returns a string reprentation of the histogram,
// which is useful for printing to a terminal.
func (h *WeightedHistogram) String() (str string) {
	str += fmt.Sprintln("Total:", h.total)

	for i := range h.bins {
		var bar string
		for j := 0; j < int(float64(h.bins[i].count)/float64(h.total)*200); j++ {
			bar += "."
		}
		str += fmt.Sprintln(h.bins[i].value, "\t", bar)
	}

	return
}

//BinsCount implements Histogram
func (h *WeightedHistogram) BinsCount() int {
	return len(h.bins)
}

//Bins implements Histogram
func (h *WeightedHistogram) Bins(i int) (float64, float64) {
	if i < 0 || i >= len(h.bins) {
		return 0, 0
	}
	b := h.bins[i]
	return b.count, b.value
}

//Alpha returns decay factor
func (h *WeightedHistogram) Alpha() float64 {
	return h.alpha
}
