extern crate num_bigint;
extern crate num_traits;
extern crate clap;
use clap::{value_t, App, Arg};
use num_bigint::{BigUint};
use std::io;
use std::io::prelude::*;
use std::collections::{VecDeque, HashSet};
use std::cmp::{min, max};

const DEFAULT_MAX_STEPS: usize = 1000;

struct Plan {
    assignment: Vec<u32>,
    width: usize,
    height: usize,
    n_districts: usize,
}

fn polyomino_splits(plan: &Plan, d1: u32, d2: u32) -> Vec<Plan> {
    let width = plan.width as usize;
    let height = plan.height as usize;
    let n_districts = plan.n_districts as usize;
    let plan_size = width * height;
    let dist_size = plan_size / n_districts;
    let mut grid = vec![false; plan_size];
    let mut dist_nodes = vec![0; 2 * dist_size]; 
    let mut node = 0;
    for (i, assign) in plan.assignment.iter().enumerate() {
        if *assign == d1 || *assign == d2 {
            grid[i] = true;
            dist_nodes[node] = i as u32;
            node += 1;
        } 
    }

    // Enumerate all possible contiguous splits of the polyomino.
    let mut seen = HashSet::<BigUint>::new();
    let mut stack = VecDeque::<Vec<u32>>::new();
    let mut districts = Vec::<Vec<u32>>::new();
    stack.push_back(vec![dist_nodes[0]]);
    seen.insert(BigUint::from(1 as u32) << dist_nodes[0]);
    while let Some(dist) = stack.pop_front() {
        let dist_bits = dist_to_bits(&dist);
        let mut neighbors = Vec::<u32>::new();
        for node in dist.iter() {
            let node_s = *node as usize;
            let y = node_s / width;
            let x = node_s % width;
            // Left neighbor
            if x > 0 && grid[node_s - 1] {
                neighbors.push(*node - 1);
            }
            // Right neighbor
            if x < width - 1 && grid[node_s + 1] {
                neighbors.push(*node + 1);
            }
            // Top neighbor
            if y > 0 && grid[node_s - width] {
                neighbors.push(*node - width as u32);
            }
            // Bottom neighbor
            if y < height - 1 && grid[node_s + width] {
                neighbors.push(*node + width as u32);
            }
        }
        for neighbor in neighbors.iter() {
            let next_bits = dist_bits.clone() | (BigUint::from(1 as u32) << neighbor);
            // TODO: uniqueness issue here?
            if !seen.contains(&next_bits) {
                let mut next_dist = dist.clone();
                next_dist.push(*neighbor);
                if dist.len() == dist_size - 1 {
                    districts.push(next_dist);
                } else {
                    // Keep building the district.
                    stack.push_back(next_dist);
                }
                seen.insert(next_bits);
            }
        }
    }

    // Check that complementing districts are contiguous.
    let mut plans = Vec::<Plan>::new();
    for district in districts.iter() {
        let mut comp = grid.clone();
        for node in district.iter() {
            comp[*node as usize] = false;
        }
        let mut start_node = dist_nodes[0];
        for node in dist_nodes.iter() {
            if comp[*node as usize] {
                start_node = *node;
                break;
            }
        }
        if comp[start_node as usize] {
            let mut comp_seen = HashSet::<u32>::new();
            let mut comp_stack = VecDeque::<u32>::new();
            comp_seen.insert(start_node);
            comp_stack.push_back(start_node);
            while let Some(node) = comp_stack.pop_front() {
                let node_s = node as usize;
                let mut neighbors = Vec::<u32>::new();
                let y = node_s / width;
                let x = node_s % width;
        
                // Left neighbor
                if x > 0 && comp[node_s - 1] {
                    neighbors.push(node - 1);
                }
                // Right neighbor
                if x < width - 1 && comp[node_s + 1] {
                    neighbors.push(node + 1);
                }
                // Top neighbor
                if y > 0 && comp[node_s - width] {
                    neighbors.push(node - width as u32);
                }
                // Bottom neighbor
                if y < height - 1 && comp[node_s + width] {
                    neighbors.push(node + width as u32);
                }
                for neighbor in neighbors.iter() {
                    if !comp_seen.contains(neighbor) {
                        comp_stack.push_back(*neighbor);
                        comp_seen.insert(*neighbor);
                    }
                }
            }
            if comp_seen.len() == dist_size {
                let mut assignment = plan.assignment.clone();
                let mut assignment_inv = plan.assignment.clone();
                for node in district.iter() {
                    assignment[*node as usize] = d1;
                    assignment_inv[*node as usize] = d2;
                }
                for node in comp_seen.into_iter() {
                    assignment[node as usize] = d2;
                    assignment_inv[node as usize] = d1;
                }

                // Make sure we didn't reproduce the current plan.
                let mut same = true;
                let mut inv_same = false;
                for node in dist_nodes.iter() {
                    if plan.assignment[*node as usize] != assignment[*node as usize] {
                        same = false;
                        break;
                    }
                }
                if !same {
                    inv_same = true;
                    for node in dist_nodes.iter() {
                        if plan.assignment[*node as usize] != assignment_inv[*node as usize] {
                            inv_same = false;
                            break;
                        }
                    }
                }

                if !same && !inv_same {
                    let next_plan = Plan {
                        assignment: assignment,
                        width: plan.width,
                        height: plan.height,
                        n_districts: plan.n_districts
                    };
                    plans.push(next_plan);
                }
            }
        }
    }
    return plans;
}

fn dist_to_bits(poly: &Vec<u32>) -> BigUint {
    let mut bits = BigUint::from(0 as u32);
    for v in poly.iter() {
        bits |= BigUint::from(1 as u32) << v;
    }
    return bits;
}

fn adjacent_dists(plan: &Plan) -> Vec<(u32, u32)> {
    let mut adj = vec![false; plan.n_districts * plan.n_districts];
    for (i, assn) in plan.assignment.iter().enumerate() {
        let y = i / plan.width;
        let x = i % plan.width;
        let mut adj_dists = vec![];
        // Left neighbor
        if x > 0 {
            adj_dists.push(plan.assignment[i - 1]);
        }
        // Right neighbor
        if x < plan.width - 1 {
            adj_dists.push(plan.assignment[i + 1]);
        }
        // Top neighbor
        if y > 0 {
            adj_dists.push(plan.assignment[i - plan.width]);
        }
        // Bottom neighbor
        if y < plan.height - 1 {
            adj_dists.push(plan.assignment[i + plan.width]);
        }
        for dist in adj_dists.iter() {
            let assn_idx = (*assn as usize) - 1;
            let dist_idx = (*dist as usize) - 1;
            adj[assn_idx * plan.n_districts + dist_idx] = true;
            adj[dist_idx * plan.n_districts + assn_idx] = true;
        }
    }
    let mut pairs = vec![];
    for i in 0..plan.n_districts {
        for j in i + 1..plan.n_districts {
            if adj[i * plan.n_districts + j] {
                pairs.push((i as u32 + 1, j as u32 + 1));
            }
        }
    }
    return pairs;
}


impl Plan {
    fn canonicalize(
        &self,
        energy_fn: &dyn Fn(&Plan) -> u64,
        freeze_fn: &dyn Fn(&Plan) -> Vec<u32>,
        max_steps: usize
    ) -> i32 {
        let mut energy = energy_fn(self);
        let mut curr = self.clone();
        let mut hist = Vec::<Plan>::new();
        let mut steps = 0;
        for _ in 0..max_steps {
            //println!("energy: {}", energy);
            if energy == 0 {
                break;
            }
            let mut splits = Vec::<Plan>::new();
            let frozen = freeze_fn(&curr);
            for (d1, d2) in adjacent_dists(&curr).iter() {
                if !frozen.contains(d1) && !frozen.contains(d2) {
                    let mut pair_splits = polyomino_splits(&curr, *d1, *d2);
                    splits.append(&mut pair_splits);
                    /*
                    for split in polyomino_splits(&curr, *d1, *d2) {
                        split.print();
                        println!();
                    }
                    */
                }
            }


            let energies: Vec<u64> = splits.iter().map(|s| energy_fn(s)).collect();
            /*
            println!("----------");
            for (i, proposal) in splits.iter().enumerate() {
                println!("energy: {}", energies[i]);
                proposal.print();
                println!();
            }
            */
            energy = *energies.iter().min().unwrap();
            for (i, split) in splits.iter().enumerate() {
                if energies[i] == energy {
                    hist.push(split.clone());
                    curr = split.clone();
                    //curr.print();
                    //println!();
                    break;
                }
            }
            steps += 1;
        }
        if energy == 0 {
            return steps;
        }
        return -1;
    }

    fn print(&self) {
        for row in 0..self.height {
            for col in 0..self.width {
                print!("{}", self.assignment[(col + (row * self.width)) as usize]);
            }
            println!();
        }
    }

    fn clone(&self) -> Plan {
        return Plan {
            assignment: self.assignment.clone(),
            width: self.width,
            height: self.height,
            n_districts: self.n_districts
        };
    }
}

fn horizontal_strip_energy(plan: &Plan) -> u64 {
    let mut horizontal_score = 0;
    for row in 0..plan.height {
        let mut consec = 1;
        for col in 1..plan.width {
            if plan.assignment[(row * plan.width) + (col - 1)] ==
               plan.assignment[(row * plan.width) + col] {
                consec += 1;
            } else {
                horizontal_score += consec * consec; //* consec;
                consec = 1;
            }
        }
        horizontal_score += consec * consec; //* consec;
    }

    /*
    let mut vertical_score = 0;
    for col in 0..plan.width {
        let mut consec = 1;
        for row in 1..plan.height {
            if plan.assignment[(col * plan.height) + (row - 1)] ==
               plan.assignment[(col * plan.height) + row] {
                consec += 1;
            } else {
                vertical_score += consec * consec;
                consec = 1;
            }
        }
        vertical_score += consec * consec;
    }
    */

    //println!("horizontal: {}\tvertical: {}", horizontal_score, vertical_score);
    //return (plan.width * plan.width * plan.width * plan.height) as u64 - horizontal_score;
    return (plan.width * plan.width * plan.height) as u64 - horizontal_score;
           //vertical_score - (plan.width * plan.height) as u64;
}

fn horizontal_strip_freeze(plan: &Plan) -> Vec<u32> {
    let mut row_min = vec![plan.height; plan.n_districts];
    let mut row_max = vec![0; plan.n_districts];
    for (i, assn) in plan.assignment.iter().enumerate() {
        let assn_idx = (*assn as usize) - 1;
        let y = i / plan.width as usize;
        row_min[assn_idx] = min(row_min[assn_idx], y);
        row_max[assn_idx] = max(row_max[assn_idx], y);
    }
    let mut frozen = vec![];
    for (i, rmin) in row_min.iter().enumerate() {
        if *rmin == row_max[i] {
            frozen.push(i as u32 + 1);
        }
    }
    return frozen;
}


fn main() {
    let matches = App::new("canonicalize")
        .version("0.1.0")
        .author("Parker J. Rule <parker.rule@tufts.edu>")
        .arg(
            Arg::with_name("width")
                .short("w")
                .long("width")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("height")
                .short("h")
                .long("height")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("n_districts")
                .short("n")
                .long("n-districts")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("max_steps")
                .short("ms")
                .long("max-steps")
                .takes_value(true),
        )
        .get_matches();
    let width = value_t!(matches.value_of("width"), usize).unwrap_or_else(|e| e.exit());
    let height = value_t!(matches.value_of("height"), usize).unwrap_or_else(|e| e.exit());
    let n_districts = value_t!(matches.value_of("n_districts"), usize).unwrap_or_else(|e| e.exit());
    let max_steps = value_t!(matches.value_of("max_steps"), usize).unwrap_or(DEFAULT_MAX_STEPS);
    //let exhaustive = matches_
    let plan_size = width * height;
    let dist_size = plan_size / n_districts;

    assert!(width > 0, "Width must be positive.");
    assert!(height > 0, "Height must be positive.");
    assert!(n_districts > 0, "Number of districts must be positive.");
    assert!(
        (width * height) % n_districts == 0,
        "Number of districts must divide grid size."
    );

    let stdin = io::stdin();
    for (i, line) in stdin.lock().lines().enumerate() {
        let plan = line.unwrap();
        let raw_assignment: Vec<String> = plan.split(',').map(|c| c.to_owned()).collect();
        let assignment: Vec<u32> = raw_assignment.iter().map(|c| c.parse().unwrap()).collect();
        if assignment.len() != plan_size {
            panic!(format!("Assignment must have length {}", plan_size));
        }
        for i in 1..n_districts + 1 {
            let dist_count = assignment.iter().filter(|&n| *n == (i as u32)).count();
            assert!(
                dist_count == dist_size, 
                format!(
                    "Expected {} units with district assignment {}",
                    dist_size, i
                )
            );
        }
        let plan = Plan {
            assignment: assignment,
            width: width,
            height: height,
            n_districts: n_districts,
        };
        let n_steps = plan.canonicalize(
            &horizontal_strip_energy,
            &horizontal_strip_freeze,
            max_steps
        );
        if n_steps > 0 {
            //println!("{} succeeded in {} steps", i, n_steps);
        } else {
            println!("{} failed", i);
        }
    }
}
