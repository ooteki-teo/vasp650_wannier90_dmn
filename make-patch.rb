#!/usr/bin/env ruby
# SPDX-License-Identifier: MIT

require 'digest'
require 'json'
require 'open3'
require 'tempfile'

if ARGV.size != 1
  warn "Usage: #{$PROGRAM_NAME} /path/to/vasp.6.5.0"
  exit 2
end

source_root = File.expand_path(ARGV.first)
recipe_root = File.join(__dir__, 'patch')
recipe = JSON.parse(File.read(File.join(recipe_root, 'vasp-6.5.0.json')))

begin
  sources = recipe.map do |relative, spec|
    original = File.binread(File.join(source_root, relative))
    source_hash = Digest::SHA256.hexdigest(original)
    if source_hash != spec.fetch('source_sha256')
      reason = source_hash == spec.fetch('patched_sha256') ? 'already patched' : 'different from the tested VASP 6.5.0 source'
      abort "#{relative}: #{reason}. Use a fresh VASP 6.5.0 source tree."
    end

    updated = original.dup
    spec.fetch('edits').reverse_each do |edit|
      insertion = edit.fetch('add', '').b
      insertion.prepend(File.binread(File.join(recipe_root, edit.fetch('include')))) if edit.key?('include')
      updated[edit.fetch('offset'), edit.fetch('delete')] = insertion
    end
    unless Digest::SHA256.hexdigest(updated) == spec.fetch('patched_sha256')
      abort "#{relative}: generated source checksum mismatch."
    end
    [relative, original, updated]
  end

  patches = sources.map do |relative, original, updated|
    Tempfile.create('dmn-before') do |before|
      Tempfile.create('dmn-after') do |after|
        before.binmode
        after.binmode
        before.write(original)
        after.write(updated)
        before.flush
        after.flush
        diff, diagnostic, status = Open3.capture3(
          'diff', '-u', '-L', "a/#{relative}", '-L', "b/#{relative}", before.path, after.path
        )
        abort "#{relative}: could not generate the diff. #{diagnostic}" unless status.exitstatus == 1
        diff
      end
    end
  end

  STDOUT.binmode
  STDOUT.write(patches.join)
rescue Errno::ENOENT => error
  warn error.message
  exit 1
end
