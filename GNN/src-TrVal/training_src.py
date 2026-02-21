from torch_geometric.loader import DataLoader 
from torch.amp import GradScaler, autocast

import torch 
import torch.optim as optim
from pathlib import Path
import numpy as np
import random

from datetime import datetime
import time
def optimized_training_loop(n_epochs, train_data, model, criterion, optimizer, device, exp, train_id, 
                          params_dir='/home/mriveraceron/glv-research/GNN-params', 
                          logs_dir = '/home/mriveraceron/glv-research/GNN-Logs',
                          seed=42, batch_size=16, use_mixed_precision=True, patience=20, 
                          save_every=50, gradient_clip_val=1.0, nlayer= 5, neurons= 10):
    
    # Section: Save lines
    train_lines = []
    now = datetime.now()
    formatted_time = now.strftime("%Y-%m-%d %H:%M:%S.%f")
    train_lines.append(
        f"{'-' * 40}\n"
        f"Starting training\n"
        f"Timestamp: {formatted_time}\n"
        f"Experiment ID: {exp}\n"
        f"Version: {train_id}\n"
        f"The number of layer used is: {nlayer}\n"
        f"The number of neurons used is: {neurons}\n"
        f"The batch size is: {batch_size}\n"
        f"The seed used is: {seed}\n"
        f"The epochs of patience is: {patience}\n"
        f"The max_norm is: {gradient_clip_val}\n"
        f"{'-' * 40}\n"
    )
    # Section: Setup
    # Set seed for reproducibility
    def set_seed(seed=42):
        # Set ALL seeds for full reproducibility
        torch.manual_seed(seed)                 # Seed CPU 
        torch.cuda.manual_seed(seed)            # Seed GPU
        np.random.seed(seed)                    # Seed numpy
        random.seed(seed)                       # Seed python random
        torch.backends.cudnn.deterministic = True   # Ensure deterministic behavior
        torch.backends.cudnn.benchmark = False      
    set_seed(seed)
    # Setup mixed precision
    scaler = GradScaler() if use_mixed_precision and device.type == 'cuda' else None
    
    # Setup DataLoader with optimizations
    def seed_worker(worker_id):
        worker_seed = torch.initial_seed() % 2**32
        np.random.seed(worker_seed)
        random.seed(worker_seed)
    
    g = torch.Generator()
    g.manual_seed(0)
    # Setup mixed precision.
    # 16-bit precision where it’s most effective and 32-bit precision where stability is crucial
    #  help(torch.utils.data.DataLoader)
    train_loader = DataLoader(
        train_data,
        batch_size=batch_size,              # Mini-batch size for training
        shuffle=True,                       # Shuffle data each epoch         
        num_workers=2,                      # Parallel data loading
        pin_memory=True,                    # Faster GPU transfer
        worker_init_fn=seed_worker,
        generator=g
    )
    
    # automatically reducing the learning rate when training gets stuck.
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, 
        mode='min',      # Monitor loss (we want it to go DOWN)
        patience=10,     # Wait 10 epochs without improvement
        factor=0.1,      # Reduce LR*factor
    )
    
    # Move model to device
    model.to(device)
    
    # Training tracking variables
    best_epoch = 1
    best_loss = float('inf')
    epochs_without_improvement = 0
    
    line1 = f"Training with batch_size={batch_size}, mixed_precision={use_mixed_precision}"
    line2 = f"DataLoader: {len(train_loader)} batches per epoch"
    train_lines.extend([line1 + '\n', line2 + '\n'])
    
    now = time.time()
    for epoch in range(1, n_epochs + 1):
        # Start epoch timer       
        start_time = time.time()
        # Training phase
        model.train()
        total_train_loss = 0.0
        batch_count = 0
        # Progress tracking
        if epoch % 10 == 1:
            print(f"Epoch {epoch}/{n_epochs}")
        
        for batch in train_loader:
            # Move mini-batch to device
            # non_blocking=True queues the transfer to happen in the background, 
            # so your program can keep running while the data is being copied.
            batch = batch.to(device, non_blocking=True)
            # Zero gradients
            optimizer.zero_grad()
            if use_mixed_precision and scaler is not None:
                #  Use 16-bytes for forward pass
                with torch.autocast(device_type=device.type):
                    out = model(batch)                          
                    loss = criterion(out, batch.y)             
                # Mixed precision backward pass
                scaler.scale(loss).backward()                   # Scale gradients UP (FP16 if  (< 6e-5) then becomes 0)
                # Gradient clipping
                # It prevents exploding gradients by capping them to a maximum value.
                # It is calculated with the norm of the gradients. which is the square root of the sum of the squares of all gradients.
                if gradient_clip_val > 0:
                    scaler.unscale_(optimizer)                  # Scale gradients back DOWN first
                    torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clip_val)   
                # Optimizer step with scaling
                scaler.step(optimizer)
                scaler.update()                   # Update scale for next iteration
            else:
                # Section: Normal precision training
                # Standard precision
                out = model(batch)
                loss = criterion(out, batch.y)
                loss.backward()
                # Gradient clipping
                if gradient_clip_val > 0:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clip_val)
                optimizer.step()
            # Accumulate loss
            total_train_loss += loss.item()
            batch_count += 1
        # Calculate average loss
        avg_train_loss = total_train_loss / batch_count
        
        # Learning rate scheduling
        scheduler.step(avg_train_loss)
        current_lr = optimizer.param_groups[0]['lr']
        
        # Section: Best model tracking
        # Track best epoch and early stopping
        if avg_train_loss < best_loss:
            best_loss = avg_train_loss
            best_epoch = epoch
            epochs_without_improvement = 0
            
            # Save best model
            best_model_path = Path(params_dir, exp) / f'{train_id}_best.pt'       # Save best model separately or overwrite
            best_model_path.parent.mkdir(parents=True, exist_ok=True)           # Create directory if not exists
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': best_loss,
                'scaler_state_dict': scaler.state_dict() if scaler else None
            }, best_model_path)
        else:
            epochs_without_improvement += 1
        
        # Section: Checkpointing
        # Calculate epoch metrics
        epoch_duration = time.time() - start_time
        
        # Logging every 10 epochs
        if epoch % 10 == 0 or epoch <= 5:  # More frequent early logging
            epoch_line = f">> Epoch {epoch:3d}, Avg Loss: {avg_train_loss:.6f}, LR: {current_lr:.2e}, Best: {best_loss:.6f}"
            duration_line = f"   Duration: {epoch_duration:.2f}s, No improvement: {epochs_without_improvement}"
            print(epoch_line)
            print(duration_line)
            train_lines.extend([epoch_line + '\n', duration_line + '\n'])
        
        # Periodic checkpointing
        if save_every > 0 and epoch % save_every == 0:
            checkpoint_path = Path(params_dir, exp) / f'{train_id}-epoch_{epoch}.pt' 
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'loss': avg_train_loss,
                'best_loss': best_loss,
                'best_epoch': best_epoch,
                'scaler_state_dict': scaler.state_dict() if scaler else None
            }, checkpoint_path)
            print(f"Checkpoint saved: {checkpoint_path.name}")
        
        # Section: Early stopping
        if patience > 0 and epochs_without_improvement >= patience:
            early_stop_line = f">> Early stopping at epoch {epoch}. No improvement for {patience} epochs."
            print(early_stop_line)
            train_lines.append(early_stop_line + '\n')
            break
    
    # Section: Final summary
    final_line = f">> Training completed. Best epoch: {best_epoch} with loss: {best_loss:.6f}"
    print(final_line)
    train_lines.append(final_line + '\n')
    
    # Save final model
    final_save_path = Path(params_dir, exp) / f'{train_id}_final.pt'
    torch.save({
        'epoch': epoch,
        'best_epoch': best_epoch,
        'best_loss': best_loss,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'scheduler_state_dict': scheduler.state_dict(),
        'final_loss': avg_train_loss,
        'scaler_state_dict': scaler.state_dict() if scaler else None
    }, final_save_path)
    
    train_duration = time.time() - now 
    duration_line = f">> Total training time: {train_duration:.2f}s"
    print(duration_line)
    train_lines.append(duration_line + '\n')

    # Section: Save training logs
    log_path = Path(logs_dir, exp) / f'{train_id}_training.log'
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, 'w') as log_file:
        log_file.writelines(train_lines)
        
    return final_save_path, best_model_path