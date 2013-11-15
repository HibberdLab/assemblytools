#!/usr/bin/env	ruby

require 'helper'

class TestRecipSnpFinder < Test::Unit::TestCase

  context "CleanFasta" do

    setup do
      @sequence1 = "ACGTACGTACGTACGT"
      @sequence2 = "ACGTACGTACGTNNNNTGCATGCATGCATGCA"
      @sequence3 = "ACGTACGTACGTNNNNGCATGCATGCATGCATNNNNTATATAGCGCGCTATATA"
      @sequence4 = "ACGTACGTACGTACGTNNNNNTGCATGCATGCATGCATGCANNNNNNN"
      @sequence5 = "ACGTACGTACGTACGTNNNN"
      @sequence6 = "NNNNNACGTACGTACGTACGTNNNN"
      @sequence7 = "ACGTACGTACNTACNTACGTACGTACGT"
    end

    should "01 return 1 sequence" do
      a = split_around_N(@sequence1,2)
      assert a.length == 1 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "02 return 2 sequence" do
      a = split_around_N(@sequence2,2)
      assert a.length == 2 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "03 return 3 sequence" do
      a = split_around_N(@sequence3,2)
      assert a.length == 3 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "04 return 2 sequence" do
      a = split_around_N(@sequence4,2)
      assert a.length == 2 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "05 return 1 sequence" do
      a = split_around_N(@sequence5,2)
      assert a.length == 1 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "06 return 1 sequence" do
      a = split_around_N(@sequence6,2)
      assert a.length == 1 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "06b return 1 sequence" do
      a = split_around_N(@sequence7,2)
      assert a.length == 1 , "length should be 1 but is #{a.length}.\n#{a}" # .each {|s| puts s}
    end

    should "07 return 1 sequence that matches" do
      assert split_around_N(@sequence1,2)[0] == "ACGTACGTACGTACGT"
    end

    should "08 return 2 sequences that matches" do
      a = split_around_N(@sequence2,2)
      t = a[0] == "ACGTACGTACGT" and a[1] == "TGCATGCATGCATGCA"
      assert  t, "#{a[0]}\n#{@sequence2.slice(0,12)}"
    end

    should "09 return 3 sequences that matches" do
      a = split_around_N(@sequence3,2)
      l = a[0] == "ACGTACGTACGT" and a[1] == "GCATGCATGCATGCAT" and a[2] == "TATATAGCGCGCTATATA"
      assert l, "#{a[0]}\n#{@sequence3.slice(0,12)}"
    end

    should "10 return 2 sequences with no Ns that match this" do
      a = split_around_N(@sequence4,2)
      l = a[0]== "ACGTACGTACGTACGT" and a[1]=="TGCATGCATGCATGCATGCA"
      assert l, "#{a[0]}\t#{a[1]}\n#{@sequence4.slice(0,16)}\t#{@sequence4.slice(21,20)}"
    end

    should "11 return 1 sequences with no Ns that match ACGTACGTACGTACGT" do
      a = split_around_N(@sequence5,2)
      l = a[0] == "ACGTACGTACGTACGT"
      assert l, "#{a[0]}\n#{@sequence5.slice(0,16)}"
    end

    should "12 return 1 sequences with no Ns that match ACGTACGTACGTACGT" do
      a = split_around_N(@sequence6,2)
      assert a[0] == "ACGTACGTACGTACGT", "there are #{a.length} sequences\n#{a}"
    end

    should "13 return 1 sequences with no Ns that match ACGTACGTACNTACNTACGTACGTACGT" do
      a = split_around_N(@sequence7,2)
      assert a[0] == "ACGTACGTACNTACNTACGTACGTACGT", "there are #{a.length} sequences\n#{a}\n[\"ACGTACGTACNTACNTACGTACGTACGT\"]"
    end
  end
end