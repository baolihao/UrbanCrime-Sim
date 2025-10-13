function neg_idx = get_neighbor_index(i, j, num_rows, num_cols, BC_type)
% function [neg_is, neg_js] = get_neighbor_index(i, j, num_rows, num_cols, type)
% set up the index for neighboring index in the order for East, South, West, and North
% for a total of 9 different groups: 4 corners, 4 sides (without corners), and 1 interior in the
% order: East side, Southeast corner, south side, southwest corner, west side, northeast corner,
% north side, northeast corner, interior

% M Zhong

% neighbors are in the order: East, South, West, North
neg_is          = [i,     i - 1, i,     i + 1];
neg_js          = [j - 1, j,     j + 1, j];
% deal with the exceptions
switch BC_type
  case 'PBC'
% due to periodicity, if the size s = (i, j) is on 4 sides or 4 corners
% the neighbor index should go around
% for i, if the index becomes 0, it should be #rows
% for i, if the index becomes > #rows, it should 1
    ind         = neg_is == 0;
    neg_is(ind) = num_rows;
    ind         = neg_is > num_rows;
    neg_is(ind) = 1;
% for j, if index becomes 0, it should be #cols
% for j, it index bcomes > #cols, it should be 1
    ind         = neg_js == 0;
    neg_js(ind) = num_cols;
    ind         = neg_js > num_cols;
    neg_js(ind) = 1;    
  case 'noFlow'
% due to periodicity, if the size s = (i, j) is on 4 sides or 4 corners
% and \partial f/\partial n = 0, so the neighbors should match
% for i, if the index becomes 0, it should be #rows
% for i, if the index becomes > #rows, it should 1
    ind         = neg_is == 0;
    neg_is(ind) = i + 1;
    ind         = neg_is > num_rows;
    neg_is(ind) = i - 1;
% for j, if index becomes 0, it should be #cols
% for j, it index bcomes > #cols, it should be 1
    ind         = neg_js == 0;
    neg_js(ind) = j + 1;
    ind         = neg_js > num_cols;
    neg_js(ind) = j - 1;
  otherwise
    error('');
end
% package them
neg_idx         = [neg_is; neg_js];
end