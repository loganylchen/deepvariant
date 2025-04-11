/*
 * Copyright 2024 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "deepvariant/channels/read_end_channel.h"
 
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"

namespace learning {
  namespace genomics {
  namespace deepvariant {
  
  void ReadEndChannel::FillReadLevelData(
      const Read& read, const DeepVariantCall& dv_call,
      const std::vector<std::string>& alt_alleles,
      std::vector<unsigned char>& read_level_data) {
    const int distance = DistanceToReadEnd(read, dv_call);
    read_level_data = std::vector<unsigned char>(
        1, ScaleColor(distance, kMaxReadEndDistance));
  }
  void ReadEndChannel::FillRefData(const std::string& ref_bases,
                                          std::vector<unsigned char>& ref_data) {
    ref_data = std::vector<unsigned char>(
        width_, static_cast<std::uint8_t>(kMaxPixelValueAsFloat));
  }
  
  // Scales an input value to pixel range 0-254.
  std::uint8_t ReadEndChannel::ScaleColor(int value, float max_val) const {
    if (static_cast<float>(value) > max_val) {
      value = max_val;
    }
    return static_cast<int>(kMaxPixelValueAsFloat *
                            ((max_val-static_cast<float>(value)) / max_val));
  }
  
 
  int ReadEndChannel::DistanceToReadEnd(const Read& read, const DeepVariantCall& dv_call) {
    const int target_pos = dv_call.variant().start();
    const int read_start = read.alignment().position().position();
    const int read_length = read.aligned_sequence().size();
    const bool is_reverse = read.alignment().position().reverse_strand();
  
    // For RNA-seq paired-end reads
    if (read.has_next_mate_position()) {
      const int mate_start = read.next_mate_position().position();
      const bool is_mate_reverse = read.next_mate_position().reverse_strand();
      const int read_5_prime_end = is_reverse ? (read_start + read_length - 1) : read_start;
      return std::abs(target_pos - read_5_prime_end);
    } else {
        return std::min(std::abs(target_pos - read_start), std::abs(target_pos -  (read_start + read_length - 1)));
    }

  }
  }  // namespace deepvariant
  }  // namespace genomics
  }  // namespace learning
