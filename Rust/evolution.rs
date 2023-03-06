extern crate rand;
extern crate rand_chacha;

use std::clone::Clone;
use std::f32::consts::PI;
use std::fmt;
use std::ops::{Index, IndexMut};
use std::process::exit;
use std::time::{Instant, Duration};
use clap::Parser;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

const EXIT_FAILURE: i32 = 1i32;

/* Structure to store data of a cell */
struct Cell {
    pos_row: f32,   // Position
    pos_col: f32,   // Position
    mov_row: f32,   // Direction of movement
    mov_col: f32,   // Direction of movement
    choose_mov: [f32; 3],   // Genes: Probabilities of 0 turning left; 1 advance; 2 turning-right
    storage: f32,   // Food/Energy stored
    age: i32,   // Number of steps that the cell has been alive
    random_seq: Box<dyn RngCore>,   // Particular random sequence generator
    alive: bool,    // Flag indicating if the cell is still alive
}

/* Make Cell struct debug-printable (with "{:?}") */
impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Cell")
         .field("alive", &self.alive)
         .field("pos_row", &self.pos_row)
         .field("pos_col", &self.pos_col)
         .field("mov_row", &self.mov_row)
         .field("mov_col", &self.mov_col)
         .field("choose_mov", &self.choose_mov)
         .field("storage", &self.storage)
         .field("age", &self.age)
         .finish()
    }
}

/* Structure for simulation statistics */
#[derive(Default)]
struct Statistics {
    history_total_cells: i32,   // Accumulated number of cells created
    history_dead_cells: i32,    // Accumulated number of dead cells
    history_max_alive_cells: i32,   // Maximum number of cells alive in a step
    history_max_new_cells: i32, // Maximum number of cells created in a step
    history_max_dead_cells: i32,    // Maximum number of cells that died in a single step
    history_max_age: i32,   // Maximum age achieved by a cell
    history_max_food: f32,  // Maximum food level in a position of the culture
}

/* Matrix structure to simplify accessing with two coordinates to a flattened array/vector */

struct Mat<T> {
   data: Vec<T>,
   rows: usize,
   cols: usize,
}

/* Immutable index operator; i.e. matrix element accesser (getter) */
impl<T> Index<(usize, usize)> for Mat<T> {
    type Output = T;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        assert!(row < self.rows, "Tried to access an out of bouds row.");
        assert!(col < self.cols, "Tried to access an out of bounds column.");
        &self.data[row * self.cols + col]
    }
}

/* Mutable index operator; i.e. matrix element setter */
impl<T> IndexMut<(usize, usize)> for Mat<T> {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        assert!(row < self.rows, "Tried to access an out of bouds row.");
        assert!(col < self.cols, "Tried to access an out of bounds column.");
        &mut self.data[row * self.cols + col]
    }
}

impl<T: Default + Clone> Mat<T> {
    fn new(rows: usize, cols: usize) -> Mat<T> {
        Mat {
            data: vec![T::default(); rows * cols],
            rows,
            cols,
        }
    }
}

/* Function: Choose a new direction of movement for a cell */
fn cell_new_direction(cell: &mut Cell) {
    let angle: f32 = 2f32 * PI * cell.random_seq.gen_range(0f32..1f32);
    cell.mov_row = angle.sin();
    cell.mov_col = angle.cos();
}

/* Function: Mutation of the movement genes on a new cell */
fn cell_mutation(cell: &mut Cell) {
    /* 1. Select which genes to change:
        0 Left grows taking part of the Advance part
        1 Advance grows taking part of the Left part
        2 Advance grows taking part of the Right part
        3 Right grows taking part of the Advance part
     */
    let mutation_type: i32 = cell.random_seq.gen_range(0i32..4i32);
    /* 2. Select the amount of mutation (up to 50%) */
    let mutation_percentge: f32 = cell.random_seq.gen_range(0f32..0.5f32);
    /* 3. Apply the mutation */
    let mutation_value: f32;
    match mutation_type{
        0 => {
            mutation_value = cell.choose_mov[1] * mutation_percentge;
            cell.choose_mov[1] -= mutation_value;
            cell.choose_mov[0] += mutation_value;
        }
        1 => {
            mutation_value = cell.choose_mov[0] * mutation_percentge;
            cell.choose_mov[0] -= mutation_value;
            cell.choose_mov[1] += mutation_value;
        }
        2 => {
            mutation_value = cell.choose_mov[2] * mutation_percentge;
            cell.choose_mov[2] -= mutation_value;
            cell.choose_mov[1] += mutation_value;
        }
        3 => {
            mutation_value = cell.choose_mov[1] * mutation_percentge;
            cell.choose_mov[1] -= mutation_value;
            cell.choose_mov[2] += mutation_value;
        }
        _ => {
            eprintln!("Error: Impossible type of mutation");
            exit(EXIT_FAILURE);
        }
    }
    /* 4. Correct potential precision problems */
    cell.choose_mov[2] = 1.0f32 - cell.choose_mov[1] - cell.choose_mov[0];
}

#[cfg(feature = "debug")]
/* Function: Print the current state of the simulation */
fn print_status(iteration: i32, culture: &Mat<f32>, num_cells: i32, cells: &Vec<Cell>, num_cells_alive: i32, sim_stat: &Statistics) {
    println!("Iteration: {}", iteration);
    print!("+");
    for _ in 0..culture.cols {
        print!("---");
    }
    println!("+");
    for i in 0..culture.rows {
        print!("|");
        for j in 0..culture.cols {
            let symbol: char;
            if culture[(i, j)] >= 20.0 {
                symbol = '+';
            } else if culture[(i, j)] >= 10.0 {
                symbol = '*';
            } else if culture[(i, j)] >= 5.0 {
                symbol = '.';
            } else {
                symbol = ' ';
            }

            let mut counter = 0i32;
            for t in 0usize..num_cells as usize {
                let row: usize = cells[t].pos_row as usize;
                let col: usize = cells[t].pos_col as usize;
                if cells[t].alive && row == i && col == j {
                    counter += 1;
                }
            }
            if counter > 9 {
                print!("(M)");
            } else if counter > 0 {
                print!("({})", counter);
            } else {
                print!(" {} ", symbol);
            }
        }
        println!("|");
    }
    print!("+");
    for _ in 0..culture.cols {
        print!("---");
    }
    println!("+");

    println!("Num_cells_alive: {:0>4}\nHistory( Cells: {:0>4}, Dead: {:0>4}, Max.alive: {:0>4}, Max.new: {:0>4}, Max.dead: {:0>4}, Max.age {:0>4}, Max.food: {:.6?}\n",
       num_cells_alive,
       sim_stat.history_total_cells,
       sim_stat.history_dead_cells,
       sim_stat.history_max_alive_cells,
       sim_stat.history_max_new_cells,
       sim_stat.history_max_dead_cells,
       sim_stat.history_max_age,
       sim_stat.history_max_food
    );
}

/* Struct: Parse the command line arguments */
#[derive(Parser,Debug)]
struct Arguments {
    /// Cultivation area size (rows)
    rows: usize,
    /// Cultivation area size (cols)
    cols: usize,
    /// Maximum number of simulation steps
    max_iter: i32,
    /// Maximum level of food on any position
    max_food: f32,
    /// Number of food sources introduced per step
    food_density: f32,
    /// Maximum number of food level in a new source
    food_level: f32,
    /// Status of the init random sequence (1)
    init_random_seq1: u16,
    /// Status of the init random sequence (2)
    init_random_seq2: u16,
    /// Status of the init random sequence (3)
    init_random_seq3: u16,
    /// Number of cells in the area at the beginning of the simulation
    num_cells: i32,
    #[clap(default_value_t=0usize)]
    /// Optional: Special food spot: Initial row
    food_spot_row: usize,
    #[clap(default_value_t=0usize)]
    /// Optional: Special food spot: Initial col
    food_spot_col: usize,
    #[clap(default_value_t=0usize)]
    /// Optional: Special food spot: Rows size
    food_spot_size_rows: usize,
    #[clap(default_value_t=0usize)]
    /// Optional: Special food spot: Cols size
    food_spot_size_cols: usize,
    #[clap(default_value_t=0f32)]
    /// Optional: Special food spot: Food density
    food_spot_denisty: f32,
    #[clap(default_value_t=0f32)]
    /// Optional: Special food spot: Food level
    food_spot_level: f32,
}

/* MAIN PROGRAM */
fn main() {
    /* 1. Read simulation arguments */
    let args = Arguments::parse();
    #[cfg(feature = "debug")]
    {
        /* Print arguments */
        println!("{:?}", args);
    }
    let food_spot_active: bool = args.food_spot_size_rows > 0 && args.food_spot_size_cols > 0; // Special food spot: Active

    let mut culture: Mat<f32> = Mat::<f32>::new(args.rows, args.cols);
    let mut culture_cells: Mat<i16> = Mat::<i16>::new(args.rows, args.cols);

    let init_random_seed: u64 = (args.init_random_seq1 as u64) << 32 | (args.init_random_seq2 as u64) << 16 | (args.init_random_seq3 as u64);
    let mut g_rng: ChaCha8Rng = ChaCha8Rng::seed_from_u64(init_random_seed);
    let mut food_random_seq: ChaCha8Rng = ChaCha8Rng::seed_from_u64(g_rng.next_u64());
    let mut food_spot_random_seq: ChaCha8Rng = ChaCha8Rng::seed_from_u64(g_rng.next_u64());

    let mut num_cells: i32 = args.num_cells;    // Number of cells currently stored in the list
    let mut cells: Vec<Cell> = Vec::<Cell>::with_capacity(num_cells as usize); // List to store cells information

    // Statistics
    let mut sim_stat: Statistics = Statistics::default();

    /* Initialize random sequences of cells */
    let cell_proto: Cell = Cell {
        pos_row: 0.0,
        pos_col: 0.0,
        mov_row: 0.0,
        mov_col: 0.0,
        choose_mov: [ 0.33f32, 0.34f32, 0.33f32 ],
        storage: 0.0,
        age: 0,
        random_seq: Box::new(g_rng.clone()),
        alive: true
    };
    for _i in 0..num_cells {
        let seed: u64 = g_rng.next_u64();
        cells.push(Cell {
            random_seq: Box::new(ChaCha8Rng::seed_from_u64(seed)),
            ..cell_proto
        });

        #[cfg(feature = "debug")]
        {
            /* Print random seed of the initial cells */
            println!("\tCell {}, Random seed: {}", _i, seed);
        }
    }

    /* 2. Start global timer */
    let timer = Instant::now();

    /* 3. Initialize culture surface and initial cells */
    for cell in cells.iter_mut() {
        // Initial age: Between 1 and 20
        cell.age = cell.random_seq.gen_range(1i32..20i32);
        // Initial storage: Between 10 and 20 units
        cell.storage = cell.random_seq.gen_range(10f32..20f32);
        // Initial position: Anywhere in the culture arena
        cell.pos_row = cell.random_seq.gen_range(0f32..args.rows as f32);
        cell.pos_col = cell.random_seq.gen_range(0f32..args.cols as f32);
        // Movement direction: Unity vector in a random direction
        cell_new_direction(cell);
    }

    // Statistics: Initialize total number of cells, and max. alive
    sim_stat.history_total_cells = num_cells;
    sim_stat.history_max_alive_cells = num_cells;

    /* Show initial cells data */
    #[cfg(feature = "debug")]
    {
        for cell in cells.iter() {
            println!("{:?}", cell);
        }
    }

    /* 4. Simulation */
    let mut current_max_food: f32 = 0f32;
    let mut num_cells_alive: i32 = num_cells;
    let mut iter: i32 = 0i32;
    for _iter in 0..args.max_iter {
        if current_max_food > args.max_food && num_cells_alive <= 0 {
            break;
        }

        let mut step_new_cells: i32 = 0i32;
        let mut step_dead_cells: i32 = 0i32;

        /* 4.1. Spreading new food */
        // Across the whole culture
        let num_new_sources: i32 = (args.rows as f32 * args.cols as f32 * args.food_density) as i32;
        for _ in 0..num_new_sources {
            let row: usize = food_random_seq.gen_range(0..args.rows) as usize;
            let col: usize = food_random_seq.gen_range(0..args.cols) as usize;
            let food: f32 = food_random_seq.gen_range(0f32..args.food_level);
            culture[(row, col)] += food;
        }
        // In the special food spot
        if food_spot_active {
            let num_new_sources: i32 = (args.food_spot_size_rows as f32 * args.food_spot_size_cols as f32 * args.food_spot_denisty) as i32;
            for _ in 0..num_new_sources {
                let row: usize = food_spot_random_seq.gen_range(args.food_spot_row..args.food_spot_size_rows) as usize;
                let col: usize = food_spot_random_seq.gen_range(args.food_spot_col..args.food_spot_size_cols) as usize;
                let food: f32 = food_spot_random_seq.gen_range(0f32..args.food_spot_level);
                culture[(row, col)] += food;
            }
        }

        /* 4.2. Prepare ancillary data structures */
        /* 4.2.1. Clear ancillary structure of the culture to account alive cells in a position after movement */
        for i in 0usize..args.rows {
            for j in 0usize..args.cols {
                culture_cells[(i, j)] = 0i16; // In the original C code, this was intentionally left as a f32
                                              // (a inoffensive "bug" for the students to find),
                                              // however, rust's compiler does not allow it.
            }
        }
        /* 4.2.2. Allocate ancillary structure to store the food level to be shared by cells in the same culture place */
        let mut food_to_share: Vec<f32> = vec![0.0f32; args.num_cells as usize];

        /* 4.3. Cells movements */
        let mut i: usize = 0usize;
        for cell in cells.iter_mut() {
            if cell.alive {
                cell.age += 1;
                // Statistics: Max age of a cell in the simulation history
                if cell.age > sim_stat.history_max_age {
                    sim_stat.history_max_age = cell.age;
                }

                /* 4.3.1. Check if the cell has the needed energy to move or stay alive */
                if cell.storage < 0.1f32 {
                    // Cell has died
                    cell.alive = false;
                    num_cells_alive -= 1;
                    step_dead_cells += 1;
                    continue;
                }
                if cell.storage < 1.0f32 {
                    // Almost dead cell, it cannot move, only if enough food is dropped here it will survive
                    cell.storage -= 0.2f32;
                } else {
                    // Consume energy to move
                    cell.storage -= 1.0f32;

                    /* 4.3.2. Choose movement direction */
                    let prob: f32 = cell.random_seq.gen_range(0f32..1f32);
                    if prob < cell.choose_mov[0] {
                        // Turn left (90 degrees)
                        let tmp: f32 = cell.mov_col;
                        cell.mov_col = cell.mov_row;
                        cell.mov_row = -tmp;
                    } else if prob >= cell.choose_mov[0] + cell.choose_mov[1] {
                        // Turn right (90 degrees)
                        let tmp: f32 = cell.mov_row;
                        cell.mov_row = cell.mov_col;
                        cell.mov_col = -tmp;
                    }
                    // else do not change direction

                    /* 4.3.3. Update position moving in the chosen direction */
                    cell.pos_row += cell.mov_row;
                    cell.pos_col += cell.mov_col;
                    // Periodic arena: Left/Right edges are connected, Top/Bottom edges are connected
                    if cell.pos_row < 0f32 { cell.pos_row += args.rows as f32; }
                    if cell.pos_row >= args.rows as f32 { cell.pos_row -= args.rows as f32; }
                    if cell.pos_col < 0f32 { cell.pos_col += args.cols as f32; }
                    if cell.pos_col >= args.cols as f32 { cell.pos_col -= args.cols as f32; }
                }

                /* 4.3.4. Annotate that there is one more cell in this culutre position */
                culture_cells[(cell.pos_row as usize, cell.pos_col as usize)] += 1;
                /* 4.3.5. Annotate the amount of food to be shared in this culture position */
                food_to_share[i] = culture[(cell.pos_row as usize, cell.pos_col as usize)];
            }
            i += 1;
        } // End of cell movements

        /* 4.4. Cell actions */
        // Space for the list of new cells (maximum number of new cells is num_cells)
        let mut new_cells: Vec<Cell> = Vec::<Cell>::with_capacity(num_cells as usize);

        i = 0usize;
        for cell in cells.iter_mut() {
            if cell.alive {
                /* 4.4.1. Food harvesting */
                let food: f32 = food_to_share[i];
                let count: i16 = culture_cells[(cell.pos_row as usize, cell.pos_col as usize)];
                let my_food: f32 = food / count as f32;
                cell.storage += my_food;

                /* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
                if cell.age > 30 && cell.storage > 20.0f32 {
                    // Split: Create new cell
                    num_cells_alive += 1;
                    sim_stat.history_total_cells += 1;
                    step_new_cells += 1;

                    new_cells.push(Cell {
                        // Random seed fpr the cell, obtained using the parent random sequence
                        random_seq: Box::new(ChaCha8Rng::seed_from_u64(cell.random_seq.next_u64())),
                        // New cell is a copy of parent cell
                        ..*cell
                    });

                    // Split energy stored and update age in both cells
                    cell.storage /= 2.0f32;
                    new_cells[step_new_cells as usize - 1].storage /= 2.0f32;
                    cell.age = 1i32;
                    new_cells[step_new_cells as usize - 1].age = 1i32;

                    // Both cells start in random directions
                    cell_new_direction(cell);
                    cell_new_direction(&mut new_cells[step_new_cells as usize - 1]);

                    // Mutations of the movement genes in both cells
                    cell_mutation(cell);
                    cell_mutation(&mut new_cells[step_new_cells as usize - 1]);
                }
            }
            i += 1;
        } // End cell actions

        /* 4.5. Clean ancillary data structures */
        /* 4.5.1. Clean the food consumed by the cells in the culture data structure */
        for cell in cells.iter() {
            if cell.alive {
                culture[(cell.pos_row as usize, cell.pos_col as usize)] = 0.0f32;
            }
        }

        /* 4.6. Clean dead cells from the original list */
        cells.retain(|cell| cell.alive);
        num_cells = cells.len() as i32;

        /* 4.7. Join cell lists: Old and new cells list */
        if step_new_cells > 0 {
            for cell in new_cells {
                cells.push(cell);
            }
        }

        /* 4.8. Decrease non-harvested food */
        current_max_food = 0.0f32;
        for i in 0..args.rows {
            for j in 0..args.cols {
                culture[(i, j)] *= 0.95f32; // Reduce 5%
                if culture[(i, j)] > current_max_food {
                    current_max_food = culture[(i, j)];
                }
            }
        }

        /* 4.9. Statistics */
        // Statistics: Max food
        if current_max_food > sim_stat.history_max_food { sim_stat.history_max_food = current_max_food; }
        // Statistics: Max new cells per step
        if step_new_cells > sim_stat.history_max_new_cells { sim_stat.history_max_new_cells = step_new_cells; }
        // Statistics: Accumulated dead and Max dead cells per step
        sim_stat.history_dead_cells += step_dead_cells;
        if step_dead_cells > sim_stat.history_max_dead_cells { sim_stat.history_max_dead_cells = step_dead_cells; }
        // Statistics: Max alive cells per step
        if num_cells_alive > sim_stat.history_max_alive_cells { sim_stat.history_max_alive_cells = num_cells_alive; }

        #[cfg(feature = "debug")]
        {
            print_status(_iter, &culture, num_cells, &cells, num_cells_alive, &sim_stat);
        }
        iter += 1;
    }

    /* 5. Stop global time */
    let ttotal: Duration = timer.elapsed();

    #[cfg(feature = "debug")]
    {
        println!("List of cells at the end of the simulation: {}", num_cells);
        for cell in cells.iter() {
            println!("{:?}", cell);
        }
    }

    /* 6. Output for leaderboard */
    println!();
    /* 6.1. Total computation time */
    println!("Time: {:?}", ttotal);

    /* 6.2. Results: Number of iterations and other statistics */
    print!("Result: {}, ", iter);
    println!("{}, {}, {}, {}, {}, {}, {} {}",
        num_cells_alive,
        sim_stat.history_total_cells,
        sim_stat.history_dead_cells,
        sim_stat.history_max_alive_cells,
        sim_stat.history_max_new_cells,
        sim_stat.history_max_dead_cells,
        sim_stat.history_max_age,
        sim_stat.history_max_food
    );

    /* 7. Free resources */

    /* 8. End */
}

