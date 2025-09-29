from torch_geometric.loader import DataLoader as GeometricDataLoader
from torch.amp import GradScaler, autocast
import time
import torch 
import torch.optim as optim
from pathlib import Path
import numpy as np
import random

def optimized_training_loop(n_epochs, train_data, model, criterion, optimizer, device, exp, train_id, 
                          save_dir='/home/mriveraceron/glv-research/GNN-params', seed=42,
                          batch_size=16, use_mixed_precision=True, patience=20, 
                          save_every=50, gradient_clip_val=1.0):
    """
    Optimized training loop with multiple improvements:
    - Batched training with DataLoader
    - Mixed precision training
    - Early stopping
    - Gradient clipping
    - Learning rate scheduling
    - Periodic checkpointing
    - Better memory management
    """
    
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
    
    # FIXME: Add seed setting for reproducibility
    # Setup DataLoader with optimizations
    def worker_init_fn(worker_id):
        np.random.seed(42 + worker_id)
    
    # Setup mixed precision.
    # 16-bit precision where it’s most effective and 32-bit precision where stability is crucial
    #  help(torch.utils.data.DataLoader)
    train_loader = GeometricDataLoader(
        train_data,
        batch_size=batch_size,
        shuffle=True,
        num_workers=2,  # Parallel data loading
        pin_memory=True,  # Faster GPU transfer
        worker_init_fn=worker_init_fn,
        generator=torch.Generator().manual_seed(42)
    )
    
    # automatically reducing the learning rate when training gets stuck.
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, 
        mode='min',      # Monitor loss (we want it to go DOWN)
        patience=10,     # Wait 10 epochs without improvement
        factor=0.5,      # Reduce LR by half (multiply by 0.5)
    )
    
    # Move model to device
    model.to(device)
    
    # Training tracking variables
    train_lines = []
    best_epoch = 1
    best_loss = float('inf')
    epochs_without_improvement = 0
    
    line1 = f"Training with batch_size={batch_size}, mixed_precision={use_mixed_precision}"
    line2 = f"DataLoader: {len(train_loader)} batches per epoch"
    train_lines.extend([line1 + '\n', line2 + '\n'])
    
    
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
            batch = batch.to(device, non_blocking=True)
            # Zero gradients
            optimizer.zero_grad()
            if use_mixed_precision and scaler is not None:
                # Mixed precision forward pass
                with autocast():
                    out = model(batch)                          #  Use 16-bytes for forward pass
                    loss = criterion(out, batch.y)              #  Use 32-bytes for loss computation
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
        
        # Track best epoch and early stopping
        if avg_train_loss < best_loss:
            best_loss = avg_train_loss
            best_epoch = epoch
            epochs_without_improvement = 0
            
            # Save best model
            best_model_path = Path(save_dir, exp) / f'{train_id}_best.pt'       # Save best model separately or overwrite
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
            checkpoint_path = Path(save_dir, exp) / f'{train_id}-epoch_{epoch}.pt' 
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
        
        # Early stopping
        if patience > 0 and epochs_without_improvement >= patience:
            early_stop_line = f">> Early stopping at epoch {epoch}. No improvement for {patience} epochs."
            print(early_stop_line)
            train_lines.append(early_stop_line + '\n')
            break
    
    # Final summary
    final_line = f">> Training completed. Best epoch: {best_epoch} with loss: {best_loss:.6f}"
    print(final_line)
    train_lines.append(final_line + '\n')
    
    # Save final model
    final_save_path = Path(save_dir) / f'{train_id}-Exp_{exp}_final.pt'
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
    
    return train_lines, final_save_path, best_model_path