import React, { useEffect, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { CiEdit } from "react-icons/ci";
import { MdDeleteOutline, MdPlayArrow, MdExpandMore, MdExpandLess, MdRefresh, MdFullscreen, MdClose, MdTextSnippet, MdTableChart, MdDelete } from "react-icons/md";
import type { IconBaseProps } from 'react-icons';
import Sidebar from '../Sidebar';
import CustomDropdown from '../CustomDropdown';
import ChemicalRenderer from '../ChemicalRenderer';
import { sessionService, workerService, renderService, type Session } from '../../services';
import './SessionDetails.css';

// Types for parsed session content
interface ParsedReaction {
  reaction_index: number;
  content: string;
  subContents: Array<Record<string, any>>;
}

interface SessionContentItem {
  type: string;
  content: string;
  aux?: {
    info_type?: string;
    reaction_index?: number;
    [key: string]: any;
  };
  timestamp?: string;
  worker_id?: string;
  final_csv?: any;
}

const SessionDetails: React.FC = () => {
  const { sessionId } = useParams<{ sessionId: string }>();
  const navigate = useNavigate();
  const [session, setSession] = useState<Session | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const [isEditing, setIsEditing] = useState(false);
  const [activeTab, setActiveTab] = useState('text');
  const [queryText, setQueryText] = useState('');
  const [isSubmittingQuery, setIsSubmittingQuery] = useState(false);
  const [selectedModel, setSelectedModel] = useState('GPT-OSS-120B');
  const [isStartingWorker, setIsStartingWorker] = useState(false);
  const [isStoppingWorker, setIsStoppingWorker] = useState(false);
  const [workerStatus, setWorkerStatus] = useState<string | null>(null);
  const [currentWorkerStatus, setCurrentWorkerStatus] = useState<any>(null);
  const [editForm, setEditForm] = useState({
    content: [] as Array<Record<string, any>>,
    state: '',
    worker_id: '',
    query: '',
    file_paths: [] as Array<Record<string, any>>
  });
  const [contentString, setContentString] = useState('');
  const [parsedContent, setParsedContent] = useState<ParsedReaction[]>([]);
  const [viewMode, setViewMode] = useState<'parsed' | 'raw'>('parsed');
  const [expandedReactions, setExpandedReactions] = useState<Set<number>>(new Set());
  const [isStreaming, setIsStreaming] = useState(false);
  const [eventSource, setEventSource] = useState<EventSource | null>(null);
  const [editableStructures, setEditableStructures] = useState<{[key: string]: string}>({});
  const [individualModal, setIndividualModal] = useState<{
    isOpen: boolean;
    structureKey: string;
    structure: string;
    type: 'SMILES' | 'SMARTS';
    index: number;
    infoType: string;
  } | null>(null);
  
  // Output files state
  const [promptFile, setPromptFile] = useState<{
    content: string | null;
    exists: boolean;
    file_path: string;
    message?: string;
  } | null>(null);
  const [csvFile, setCsvFile] = useState<{
    data: Array<Record<string, any>> | null;
    headers: string[] | null;
    exists: boolean;
    file_path: string;
    row_count: number;
    message?: string;
  } | null>(null);
  const [loadingOutputFiles, setLoadingOutputFiles] = useState(false);
  
  // Modal states for output files
  const [promptModal, setPromptModal] = useState<{
    isOpen: boolean;
    content: string;
  }>({ isOpen: false, content: '' });
  const [csvModal, setCsvModal] = useState<{
    isOpen: boolean;
    data: Array<Record<string, any>>;
    headers: string[];
    isEditing: boolean;
  }>({ isOpen: false, data: [], headers: [], isEditing: false });
  
  // CSV editing state
  const [editedCsvData, setEditedCsvData] = useState<Array<Record<string, any>>>([]);
  const [isSavingCsv, setIsSavingCsv] = useState(false);
  const [isMergingCsv, setIsMergingCsv] = useState(false);

  // Modal control functions
  const openIndividualModal = (structureKey: string, structure: string, type: 'SMILES' | 'SMARTS', index: number, infoType: string) => {
    setIndividualModal({ isOpen: true, structureKey, structure, type, index, infoType });
  };

  const closeIndividualModal = () => {
    setIndividualModal(null);
  };

  // Output file modal controls
  const openPromptModal = (content: string) => {
    setPromptModal({ isOpen: true, content });
  };

  const closePromptModal = () => {
    setPromptModal({ isOpen: false, content: '' });
  };

  const openCsvModal = (data: Array<Record<string, any>>, headers: string[]) => {
    setCsvModal({ isOpen: true, data, headers, isEditing: false });
    setEditedCsvData(JSON.parse(JSON.stringify(data))); // Deep copy for editing
  };

  const closeCsvModal = () => {
    setCsvModal({ isOpen: false, data: [], headers: [], isEditing: false });
    setEditedCsvData([]);
    setIsSavingCsv(false);
  };

  const toggleCsvEditing = () => {
    setCsvModal(prev => ({ ...prev, isEditing: !prev.isEditing }));
  };

  const handleCsvCellChange = (rowIndex: number, column: string, value: string) => {
    setEditedCsvData(prev => {
      const newData = [...prev];
      newData[rowIndex] = { ...newData[rowIndex], [column]: value };
      return newData;
    });
  };

  const deleteCsvRow = (rowIndex: number) => {
    setEditedCsvData(prev => {
      const newData = [...prev];
      newData.splice(rowIndex, 1);
      return newData;
    });
  };

  const saveCsvChanges = async () => {
    if (!sessionId) return;
    
    try {
      setIsSavingCsv(true);
      
      // Save the edited data to the backend
      const result = await sessionService.saveSessionCsvData(sessionId, editedCsvData);
      
      // Update the local state with the saved data
      setCsvModal(prev => ({ ...prev, data: editedCsvData, isEditing: false }));
      setCsvFile(prev => prev ? { 
        ...prev, 
        data: editedCsvData,
        row_count: editedCsvData.length 
      } : null);
      
      console.log('CSV data saved successfully:', result);
      
    } catch (error) {
      console.error('Failed to save CSV data:', error);
      // Keep edit mode active on error
      alert('Failed to save CSV data. Please try again.');
    } finally {
      setIsSavingCsv(false);
    }
  };

  const resetCsvChanges = () => {
    setEditedCsvData(JSON.parse(JSON.stringify(csvModal.data)));
  };

  const mergeCsvToProject = async () => {
    if (!sessionId) return;
    
    try {
      setIsMergingCsv(true);
      
      // Merge the session CSV with the project master CSV
      const result = await sessionService.mergeSessionCsvToProject(sessionId);
      
      console.log('CSV merge completed successfully:', result);
      
      // Show success message with merge details
      const mergeResults = result.merge_results;
      alert(
        `CSV merge completed successfully!\n\n` +
        `Rows read from session: ${mergeResults.rows_read}\n` +
        `Rows added to project: ${mergeResults.rows_added}\n` +
        `Rows skipped (duplicates): ${mergeResults.rows_skipped}\n` +
        `Total rows in project CSV: ${mergeResults.total_rows_after_merge}`
      );
      
    } catch (error) {
      console.error('Failed to merge CSV files:', error);
      alert('Failed to merge CSV files. Please try again.');
    } finally {
      setIsMergingCsv(false);
    }
  };

  // Handle structure editing
  const handleStructureChange = (structureKey: string, newValue: string) => {
    setEditableStructures(prev => ({
      ...prev,
      [structureKey]: newValue
    }));
  };

  // Get the current value for a structure (either edited or original)
  const getCurrentStructureValue = (originalStructure: string, index: number, infoType: string) => {
    const structureKey = `${infoType}-${index}`;
    return editableStructures[structureKey] ?? originalStructure;
  };

  // Fetch output files for completed sessions
  const fetchOutputFiles = async (sessionId: string) => {
    if (!sessionId) return;
    
    try {
      setLoadingOutputFiles(true);
      
      // Fetch both prompt file and CSV file in parallel
      const [promptResponse, csvResponse] = await Promise.allSettled([
        sessionService.getSessionPromptFile(sessionId),
        sessionService.getSessionCsvFile(sessionId)
      ]);
      
      // Handle prompt file response
      if (promptResponse.status === 'fulfilled') {
        setPromptFile({
          content: promptResponse.value.content,
          exists: promptResponse.value.exists,
          file_path: promptResponse.value.file_path,
          message: promptResponse.value.message
        });
      } else {
        console.warn('Failed to fetch prompt file:', promptResponse.reason);
        setPromptFile(null);
      }
      
      // Handle CSV file response
      if (csvResponse.status === 'fulfilled') {
        setCsvFile({
          data: csvResponse.value.data,
          headers: csvResponse.value.headers,
          exists: csvResponse.value.exists,
          file_path: csvResponse.value.file_path,
          row_count: csvResponse.value.row_count,
          message: csvResponse.value.message
        });
      } else {
        console.warn('Failed to fetch CSV file:', csvResponse.reason);
        setCsvFile(null);
      }
      
    } catch (error) {
      console.error('Error fetching output files:', error);
    } finally {
      setLoadingOutputFiles(false);
    }
  };

  // Parser function to transform session content
  // Helper function to detect and extract chemical data
  const getChemicalData = (item: any) => {
    if (!item || typeof item !== 'object') return null;
    
    // Check if the item has aux.info_type indicating chemical data
    if (item.aux?.info_type) {
      const infoType = item.aux.info_type;
      
      // Extract values array from different sources
      let values: string[] = [];
      try {
        // First, check if values are in aux.returnable.values
        if (item.aux?.returnable?.values && Array.isArray(item.aux.returnable.values)) {
          values = item.aux.returnable.values;
        }
        // Then check content field
        else if (typeof item.content === 'string') {
          const parsed = JSON.parse(item.content);
          // Handle different JSON structures
          if (Array.isArray(parsed.values)) {
            values = parsed.values;
          } else if (Array.isArray(parsed)) {
            // If the parsed content is directly an array
            values = parsed;
          } else if (parsed.smiles && Array.isArray(parsed.smiles)) {
            // Handle case where values are in a smiles key
            values = parsed.smiles;
          } else if (parsed.smarts && Array.isArray(parsed.smarts)) {
            // Handle case where values are in a smarts key
            values = parsed.smarts;
          }
        } else if (item.content?.values && Array.isArray(item.content.values)) {
          values = item.content.values;
        } else if (Array.isArray(item.content)) {
          // If content is directly an array
          values = item.content;
        }
        
        // Filter out empty or invalid values
        values = values.filter(value => value && typeof value === 'string' && value.trim().length > 0);
        
      } catch (error) {
        console.warn('Failed to parse chemical content:', error);
        return null;
      }
      
      // Map info_type to render type
      if (infoType === 'smiles_update') {
        return {
          type: 'SMILES' as const,
          values: values,
          infoType: infoType
        };
      } else if (infoType === 'smarts_rxn_update' || infoType === 'smirks_rxn_update') {
        return {
          type: 'SMARTS' as const,
          values: values,
          infoType: infoType
        };
      }
    }
    
    return null;
  };

  // Helper function to render chemical structures
  const renderChemicalStructures = (chemicalData: {type: 'SMILES' | 'SMARTS', values: string[], infoType: string}) => {
    // Map info_type to display header
    const getHeaderTitle = (infoType: string) => {
      switch (infoType) {
        case 'smiles_update':
          return 'Composition SMILES';
        case 'smarts_rxn_update':
          return 'Reaction SMARTS';
        case 'smirks_rxn_update':
          return 'Reaction SMIRKS';
        default:
          return infoType.replace(/_/g, ' ').toUpperCase();
      }
    };

    if (!chemicalData.values || chemicalData.values.length === 0) {
      return (
        <div className="chemical-structures">
          <div className="chemical-header">
            <div className="chemical-header-left">
              <h4>{getHeaderTitle(chemicalData.infoType)}</h4>
              <span className="chemical-count">No structures</span>
            </div>
          </div>
          <p className="no-structures">No chemical structures found</p>
        </div>
      );
    }
    
    return (
      <div className="chemical-structures">
        <div className="chemical-header">
          <div className="chemical-header-left">
            <h4>{getHeaderTitle(chemicalData.infoType)}</h4>
            <span className="chemical-count">{chemicalData.values.length} structure{chemicalData.values.length !== 1 ? 's' : ''}</span>
          </div>
        </div>
        <div className="chemical-grid">
          {chemicalData.values.map((structure, index) => {
            const structureKey = `${chemicalData.infoType}-${index}`;
            const currentValue = getCurrentStructureValue(structure, index, chemicalData.infoType);
            
            return (
              <div key={`${chemicalData.infoType}-${index}`} className="chemical-item">
                <div className="chemical-item-header">
                  <span className="chemical-index">#{index + 1}</span>
                  <button 
                    className="individual-fullscreen-btn"
                    onClick={() => openIndividualModal(structureKey, currentValue, chemicalData.type, index, chemicalData.infoType)}
                    title="Open in fullscreen"
                  >
                    {React.createElement(MdFullscreen as React.ComponentType<IconBaseProps>, { size: 14 })}
                  </button>
                </div>
                <div className="chemical-structure-input-container">
                  <textarea
                    className="chemical-structure-textarea"
                    value={currentValue}
                    onChange={(e) => handleStructureChange(structureKey, e.target.value)}
                    placeholder={`Enter ${chemicalData.type} structure`}
                    title="Edit this chemical structure"
                    rows={2}
                  />
                </div>
                <ChemicalRenderer
                  data={currentValue}
                  inputType={chemicalData.type}
                  className="inline"
                  onRenderError={(error) => console.warn(`Failed to render ${chemicalData.type}:`, currentValue, error)}
                />
              </div>
            );
          })}
        </div>
      </div>
    );
  };

  const parseSessionContent = (content: SessionContentItem[]): ParsedReaction[] => {
    const reactions: ParsedReaction[] = [];
    
    if (!content || !Array.isArray(content)) {
      console.log('SessionDetails: No content or content is not an array');
      return reactions;
    }

    console.log('SessionDetails: Parsing session content:', content.length, 'items');

    content.forEach((item, itemIndex) => {
      if (!item || typeof item !== 'object') {
        console.log(`SessionDetails: Item ${itemIndex} is not a valid object:`, item);
        return;
      }

      console.log(`SessionDetails: Processing item ${itemIndex}:`, item);

      // Check if item has "aux" key
      if (item.aux) {
        // Handle reaction_text type - this creates a new reaction
        if (item.aux.info_type === 'reaction_text' && typeof item.aux.reaction_index === 'number') {
          const reactionIndex = item.aux.reaction_index;
          console.log(`SessionDetails: Found reaction_text for reaction ${reactionIndex}:`, item.content);
          
          // Check if this reaction already exists
          const existingReaction = reactions.find(r => r.reaction_index === reactionIndex);
          if (!existingReaction) {
            reactions.push({
              reaction_index: reactionIndex,
              content: item.content || '',
              subContents: []
            });
            console.log(`SessionDetails: Created new reaction ${reactionIndex}`);
          }
        }
        
        // Handle all other items that have reaction_index in aux - add them to subContents
        else if (typeof item.aux.reaction_index === 'number') {
          const reactionIndex = item.aux.reaction_index;
          const reaction = reactions.find(r => r.reaction_index === reactionIndex);
          
          console.log(`SessionDetails: Adding subContent to reaction ${reactionIndex}:`, item.aux.info_type);
          
          // Skip smiles_rxn_update items entirely
          if (item.aux.info_type === 'smiles_rxn_update') {
            console.log(`SessionDetails: Skipping smiles_rxn_update item for reaction ${reactionIndex}`);
            return;
          }
          
          // Only add items with type "data" to subContents
          if (item.type === 'data') {
            if (reaction) {
              reaction.subContents.push(item);
            } else {
              // If reaction doesn't exist yet, create it with empty content
              // This handles cases where subContent comes before reaction_text
              reactions.push({
                reaction_index: reactionIndex,
                content: `Reaction ${reactionIndex}`, // Default content
                subContents: [item]
              });
              console.log(`SessionDetails: Created new reaction ${reactionIndex} from subContent`);
            }
          }
        }
      }
      
      // Handle items without aux or reaction_index (like run_complete)
      else {
        console.log(`SessionDetails: Adding item without aux to general section:`, item);
        
        // Only add data type items
        if (item.type === 'data') {
          // Add to the last reaction if it exists, otherwise create a general reaction
          if (reactions.length > 0) {
            reactions[reactions.length - 1].subContents.push(item);
          } else {
            // Create a general reaction for items without reaction_index
            reactions.push({
              reaction_index: -1, // Use -1 for general items
              content: 'General Session Items',
              subContents: [item]
            });
          }
        }
      }
    });

    // Sort reactions by reaction_index (put general items at the end)
    reactions.sort((a, b) => {
      if (a.reaction_index === -1) return 1;
      if (b.reaction_index === -1) return -1;
      return a.reaction_index - b.reaction_index;
    });

    console.log('SessionDetails: Final parsed reactions:', reactions);
    return reactions;
  };

  const toggleReactionExpansion = (reactionIndex: number) => {
    setExpandedReactions(prev => {
      const newSet = new Set(prev);
      if (newSet.has(reactionIndex)) {
        newSet.delete(reactionIndex);
      } else {
        newSet.add(reactionIndex);
      }
      return newSet;
    });
  };

  const toggleAllReactions = () => {
    if (expandedReactions.size === parsedContent.length) {
      // All are expanded, collapse all
      setExpandedReactions(new Set());
    } else {
      // Some or none are expanded, expand all
      const allReactionIndexes = new Set(parsedContent.map(r => r.reaction_index));
      setExpandedReactions(allReactionIndexes);
    }
  };

  useEffect(() => {
    if (sessionId) loadSession(sessionId);
  }, [sessionId]);

  // Cleanup EventSource on unmount
  useEffect(() => {
    return () => {
      if (eventSource) {
        eventSource.close();
      }
    };
  }, [eventSource]);

  const loadSession = async (id: string) => {
    try {
      setLoading(true);
      setError('');
      const data = await sessionService.getSession(id);
      setSession(data);
      const contentStr = JSON.stringify(data.content || [], null, 2);
      setContentString(contentStr);
      setQueryText(data.query || '');
      setEditForm({
        content: data.content || [],
        state: data.state || '',
        worker_id: data.worker_id,
        query: data.query || '',
        file_paths: data.file_paths || []
      });
      
      // Parse the session content
      const parsed = parseSessionContent((data.content || []) as SessionContentItem[]);
      setParsedContent(parsed);
      
      // Auto-expand all reactions when content is first loaded
      if (parsed.length > 0) {
        const allReactionIndexes = new Set(parsed.map(r => r.reaction_index));
        setExpandedReactions(allReactionIndexes);
      }
      
      // Fetch output files if session has file_paths or is completed
      if ((data.file_paths && data.file_paths.length > 0) || 
          data.state === 'completed' || 
          data.state === 'finished' ||
          (data.content && data.content.length > 0)) {
        await fetchOutputFiles(id);
      }
    } catch (err) {
      setError('Session not found');
    } finally {
      setLoading(false);
    }
  };

  const handleMenuSelect = (menu: string) => {
    if (menu === 'projects') navigate('/');
    else if (menu === 'dashboard') navigate('/?section=dashboard');
    else if (menu === 'settings') navigate('/?section=settings');
  };
  const handleSidebarToggle = (collapsed: boolean) => setSidebarCollapsed(collapsed);

  const handleEdit = () => setIsEditing(true);
  const handleCancelEdit = () => {
    setIsEditing(false);
    if (session) {
      const contentStr = JSON.stringify(session.content || [], null, 2);
      setContentString(contentStr);
      setEditForm({ 
        content: session.content || [], 
        state: session.state || '', 
        worker_id: session.worker_id,
        query: session.query || '',
        file_paths: session.file_paths || []
      });
    }
  };
  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement | HTMLSelectElement>) => {
    const { name, value } = e.target;
    if (name === 'content') {
      setContentString(value);
      try {
        const parsedContent = JSON.parse(value);
        setEditForm(prev => ({ ...prev, content: parsedContent }));
      } catch {
        // Invalid JSON, keep the current content in editForm
      }
    } else {
      setEditForm(prev => ({ ...prev, [name]: value }));
    }
  };
  const handleSaveEdit = async () => {
    if (!session || !sessionId) return;
    try {
      setLoading(true);
      const updated = await sessionService.updateSession(sessionId, editForm);
      setSession(updated);
      setIsEditing(false);
      setError('');
    } catch {
      setError('Failed to update session');
    } finally {
      setLoading(false);
    }
  };
  const handleDelete = async () => {
    if (!session || !sessionId) return;
    if (!window.confirm('Delete this session?')) return;
    try {
      setLoading(true);
      await sessionService.deleteSession(sessionId);
      navigate(-1);
    } catch {
      setError('Failed to delete session');
      setLoading(false);
    }
  };

  const handleQuerySubmit = async () => {
    if (!sessionId || !queryText.trim()) return;
    try {
      setIsSubmittingQuery(true);
      setError('');
      await sessionService.updateSessionQuery(sessionId, queryText.trim());
      // Reload session to get updated data
      await loadSession(sessionId);
    } catch (err) {
      setError('Failed to update query');
    } finally {
      setIsSubmittingQuery(false);
    }
  };

  const startStreaming = async (sessionId: string) => {
    if (eventSource) {
      eventSource.close();
    }
    
    setIsStreaming(true);
    const clientId = `client_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
    
    console.log('🚀 Starting streaming for session:', sessionId, 'with client:', clientId);
    
    try {
      const newEventSource = await sessionService.streamSessionContent(
        sessionId,
        clientId,
        (data) => {
          // Handle incoming streaming data
          console.log(`📥 Received streaming data:`, data);
          
          if (data.type === 'data' && data.content) {
            // Update session content incrementally
            setSession(prevSession => {
              if (!prevSession) return prevSession;
              
              const currentContent = prevSession.content || [];
              const updatedContent = [...currentContent, data.content];
              
              // Update content string for raw view
              setContentString(JSON.stringify(updatedContent, null, 2));
              
              // Parse and update parsed content
              const parsed = parseSessionContent(updatedContent as SessionContentItem[]);
              setParsedContent(parsed);
              
              return {
                ...prevSession,
                content: updatedContent
              };
            });
          } else if (data.type === 'info') {
            console.log('Stream info:', data.message);
          }
        },
        (error) => {
          console.error('❌ Streaming error:', error);
          setError(`Streaming error: ${error.message}`);
          setIsStreaming(false);
        },
        () => {
          console.log('✅ Streaming completed');
          setIsStreaming(false);
          // Reload session to get final state
          loadSession(sessionId);
        }
      );
      
      setEventSource(newEventSource);
    } catch (error) {
      console.error('Failed to start streaming:', error);
      setError('Failed to start content streaming');
      setIsStreaming(false);
    }
  };

  const handleStartWorker = async () => {
    if (!sessionId || !selectedModel) return;
    try {
      setIsStartingWorker(true);
      setError('');
      
      // Only use the saved session query, not the textarea content
      const queryToUse = session?.query?.trim();
      
      const response = await workerService.startSessionWorker(sessionId, {
        model_name: selectedModel,
        query: queryToUse || undefined
      });
      
      setWorkerStatus('starting');
      // Start streaming content after worker is started
      await startStreaming(sessionId);
      // Reload session to get updated state
      await loadSession(sessionId);
    } catch (err) {
      setError('Failed to start worker');
    } finally {
      setIsStartingWorker(false);
    }
  };

  const handleRerunSession = async () => {
    if (!sessionId || !selectedModel) return;
    try {
      setIsStartingWorker(true);
      setError('');
      
      // Clear existing content before rerunning
      setSession(prevSession => {
        if (!prevSession) return prevSession;
        return {
          ...prevSession,
          content: []
        };
      });
      setContentString('[]');
      setParsedContent([]);
      
      // Close existing stream if any
      if (eventSource) {
        eventSource.close();
        setEventSource(null);
      }
      
      // Only use the saved session query
      const queryToUse = session?.query?.trim();
      
      const response = await workerService.startSessionWorker(sessionId, {
        model_name: selectedModel,
        query: queryToUse || undefined
      });
      
      setWorkerStatus('starting');
      // Start streaming content after worker is started
      await startStreaming(sessionId);
    } catch (err) {
      setError('Failed to rerun session');
    } finally {
      setIsStartingWorker(false);
    }
  };

  const handleStopWorker = async () => {
    if (!sessionId) return;
    
    console.log('🛑 Attempting to stop worker for session:', sessionId);
    
    try {
      setIsStoppingWorker(true);
      setError('');
      
      // First, call the backend to stop the worker
      console.log('📡 Calling backend to stop worker...');
      const result = await sessionService.stopSessionWorker(sessionId);
      console.log('✅ Backend stop response:', result);
      
      // Then close the streaming connection
      if (eventSource) {
        console.log('🔌 Closing EventSource connection...');
        eventSource.close();
        setEventSource(null);
      }
      
      // Stop streaming and reset states
      setIsStreaming(false);
      setIsStartingWorker(false);
      setWorkerStatus('stopped');
      
      // Reload session to get updated state
      console.log('🔄 Reloading session state...');
      await loadSession(sessionId);
      console.log('✅ Worker stopped successfully');
    } catch (err: any) {
      const errorMessage = err?.message || 'Failed to stop worker';
      setError(errorMessage);
      console.error('❌ Error stopping worker:', err);
    } finally {
      setIsStoppingWorker(false);
    }
  };

  return (
    <div className="session-details-page">
      <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu="projects" />
      <div className={`session-details-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
        <div className="breadcrumb">
          <span onClick={() => navigate('/')} className="breadcrumb-link">Projects</span>
          <span className="breadcrumb-separator">/</span>
          <span className="breadcrumb-link" onClick={() => navigate(-1)}>Project</span>
          <span className="breadcrumb-separator">/</span>
          <span className="breadcrumb-current">{session?.name || sessionId}</span>
        </div>
        {loading ? (
          <div className="session-details-loading"><div className="loading-spinner"></div><p>Loading session...</p></div>
        ) : error ? (
          <div className="session-details-error"><h2>Error</h2><p>{error}</p></div>
        ) : !session ? (
          <div className="session-details-error"><h2>Not Found</h2><p>Session not found.</p></div>
        ) : (
          <>
            <div className="session-header-bar">
              <div className="session-name">
                <h1>{session.name}</h1>
              </div>
              <div className="session-state-tag">
                <span className="state-badge">{session.state || 'No state'}</span>
              </div>
              <div className="session-actions">
                {!isEditing && (
                  <>
                    <button 
                      onClick={handleEdit} 
                      className="btn-icon btn-primary" 
                      data-tooltip="Edit Session"
                    >
                      {React.createElement(CiEdit as React.ComponentType<IconBaseProps>, { size: 20 })}
                    </button>
                    <button 
                      onClick={handleDelete} 
                      className="btn-icon btn-danger" 
                      data-tooltip="Delete Session"
                    >
                      {React.createElement(MdDeleteOutline as React.ComponentType<IconBaseProps>, { size: 20 })}
                    </button>
                  </>
                )}
              </div>
            </div>
            
            <div className="query-section">
              <div className="query-tabs">
                <button 
                  className={`query-tab ${activeTab === 'text' ? 'active' : ''}`}
                  onClick={() => setActiveTab('text')}
                >
                  Text
                </button>
                <button 
                  className={`query-tab ${activeTab === 'url' ? 'active' : ''} disabled`}
                  onClick={() => setActiveTab('url')}
                  disabled
                >
                  URL
                </button>
                <button 
                  className={`query-tab ${activeTab === 'pdf' ? 'active' : ''} disabled`}
                  onClick={() => setActiveTab('pdf')}
                  disabled
                >
                  PDF
                </button>
              </div>
              
              <div className="query-content">
                {activeTab === 'text' && (
                  <div className="text-query-panel">
                    <textarea 
                      className="query-textarea"
                      placeholder="Enter your text query here..."
                      rows={6}
                      value={queryText}
                      onChange={(e) => setQueryText(e.target.value)}
                    />
                    <div className="query-actions">
                      <button 
                        className="btn-primary query-submit-btn"
                        onClick={handleQuerySubmit}
                        disabled={isSubmittingQuery || !queryText.trim()}
                      >
                        {isSubmittingQuery 
                          ? 'Processing...' 
                          : (session?.query && session.query.trim() !== '') 
                            ? 'Update Query' 
                            : 'Submit Query'
                        }
                      </button>
                    </div>
                  </div>
                )}
                {activeTab === 'url' && (
                  <div className="url-query-panel">
                    <p className="coming-soon">URL query functionality coming soon...</p>
                  </div>
                )}
                {activeTab === 'pdf' && (
                  <div className="pdf-query-panel">
                    <p className="coming-soon">PDF query functionality coming soon...</p>
                  </div>
                )}
              </div>
            </div>
            
            <div className="runner-section">
              <div className="runner-content">
                {/* Show runner content for all states, not just init */}
                <div className="runner-init">
                  {/* Show which query will be used */}
                  <div className="query-preview">
                    <label>Query to process:</label>
                    {(() => {
                      // Only use the saved session query, not the textarea content
                      const currentQuery = session?.query?.trim();
                      if (!currentQuery) {
                        return (
                          <>
                            <p className="query-text no-query">No query saved</p>
                            <small className="query-warning">⚠️ No query saved. Please submit a query above before starting the runner.</small>
                          </>
                        );
                      }
                      
                      const maxLength = 100;
                      const isLong = currentQuery.length > maxLength;
                      
                      let displayQuery;
                      if (isLong) {
                        // Show start and end with "..." in between
                        const startLength = Math.floor((maxLength - 3) / 2); // Reserve 3 chars for "..."
                        const endLength = maxLength - 3 - startLength;
                        const start = currentQuery.substring(0, startLength);
                        const end = currentQuery.substring(currentQuery.length - endLength);
                        displayQuery = `${start}...${end}`;
                      } else {
                        displayQuery = currentQuery;
                      }
                      
                      return (
                        <>
                          <p className="query-text">
                            {displayQuery}
                          </p>
                          <small className="query-info">
                            {currentQuery.length} characters
                          </small>
                        </>
                      );
                    })()}
                  </div>

                  <div>
                    <CustomDropdown
                      label="Select Model:"
                      options={workerService.getAvailableModels()}
                      value={selectedModel}
                      onChange={setSelectedModel}
                      disabled={isStartingWorker || isStreaming} // Disable when running
                    />
                    
                    {/* Show different button states based on running state */}
                    {isStartingWorker || isStreaming ? (
                      // Running state - show Stop Runner button
                      <button 
                        className="start-runner-btn stop-runner-btn"
                        onClick={handleStopWorker}
                        disabled={isStoppingWorker}
                      >
                        {React.createElement(MdClose as React.ComponentType<IconBaseProps>, { size: 18 })}
                        <span>{isStoppingWorker ? 'Stopping...' : 'Stop Runner'}</span>
                        <span className="runner-label">Running</span>
                      </button>
                    ) : (session?.content && session.content.length > 0) || session?.state === 'cancelled' ? (
                      // Completed or Cancelled state - show Re-run button
                      <button 
                        className={`start-runner-btn rerun-runner-btn ${session?.state === 'cancelled' ? 'cancelled-state' : ''}`}
                        onClick={handleRerunSession}
                        disabled={!selectedModel || !session?.query?.trim()}
                      >
                        {React.createElement(MdRefresh as React.ComponentType<IconBaseProps>, { size: 18 })}
                        <span>Re-run</span>
                        <span className="runner-label">{session?.state === 'cancelled' ? 'Cancelled' : 'Completed'}</span>
                      </button>
                    ) : (
                      // Init state - show Start Runner button
                      <button 
                        className="start-runner-btn"
                        onClick={handleStartWorker}
                        disabled={!selectedModel || !session?.query?.trim()}
                      >
                        {React.createElement(MdPlayArrow as React.ComponentType<IconBaseProps>, { size: 18 })}
                        <span>Start Runner</span>
                        <span className="runner-label">Init</span>
                      </button>
                    )}
                  </div>
                </div>
              </div>
            </div>
            
            {/* Session Content Display */}
            <div className="content-section">
              <div className="content-header">
                <h2>Session Content</h2>
                <div className="content-header-right">
                  <div className="view-mode-toggle">
                    <button 
                      className={`view-mode-btn ${viewMode === 'parsed' ? 'active' : ''}`}
                      onClick={() => setViewMode('parsed')}
                    >
                      Parsed View
                    </button>
                    <button 
                      className={`view-mode-btn ${viewMode === 'raw' ? 'active' : ''}`}
                      onClick={() => setViewMode('raw')}
                    >
                      Raw JSON
                    </button>
                  </div>
                  {currentWorkerStatus && (
                    <div className="worker-status-indicator">
                      <span className={`status-dot ${currentWorkerStatus.status}`}></span>
                      <span className="status-text">
                        Worker: {currentWorkerStatus.status}
                      </span>
                    </div>
                  )}
                  {isStreaming && (
                    <div className="streaming-indicator">
                      <span className="streaming-dot"></span>
                      <span className="streaming-text">Live streaming...</span>
                    </div>
                  )}
                  <div className="reactions-count-indicator">
                    <span className="reactions-count-label">Reactions Found:</span>
                    <span className="reactions-count-value">
                      {parsedContent.filter(r => r.reaction_index !== -1).length}
                    </span>
                  </div>
                </div>
              </div>
              
              <div className="content-display">
                {session?.content && session.content.length > 0 ? (
                  <div className="content-container">
                    {viewMode === 'parsed' ? (
                      <div className="parsed-content">
                        {parsedContent.length > 0 ? (
                          <div className="reactions-container">
                            {parsedContent.map((reaction, index) => {
                              const isExpanded = expandedReactions.has(reaction.reaction_index);
                              return (
                              <div key={`reaction-${reaction.reaction_index}-${index}`} className="reaction-item">
                                <div className={`reaction-header ${isExpanded ? 'expanded' : ''}`} onClick={() => toggleReactionExpansion(reaction.reaction_index)}>
                                  <div className="reaction-header-content">
                                    <h3>
                                      {reaction.reaction_index === -1 
                                        ? 'General Session Information' 
                                        : `Reaction ${reaction.reaction_index}`
                                      }
                                    </h3>
                                    <div className="reaction-content-text">
                                      {reaction.content}
                                    </div>
                                  </div>
                                  <div className="reaction-expand-icon">
                                    {isExpanded ? (
                                      React.createElement(MdExpandLess as React.ComponentType<IconBaseProps>, { size: 24 })
                                    ) : (
                                      React.createElement(MdExpandMore as React.ComponentType<IconBaseProps>, { size: 24 })
                                    )}
                                  </div>
                                </div>
                                
                                {isExpanded && reaction.subContents.length > 0 && (
                                  <div className="sub-contents">
                                    <h4>
                                      {reaction.reaction_index === -1 
                                        ? 'Session Information:' 
                                        : 'Data Found:'
                                      }
                                    </h4>
                                    <div className="sub-contents-list">
                                      {reaction.subContents.map((subContent, subIndex) => (
                                        <div key={`sub-${reaction.reaction_index}-${subIndex}`} className="sub-content-item">
                                          <div className="sub-content-body">
                                            {(() => {
                                              const chemicalData = getChemicalData(subContent);
                                              
                                              // Debug logging
                                              console.log('SubContent processing:', {
                                                type: subContent.type,
                                                infoType: subContent.aux?.info_type,
                                                content: subContent.content,
                                                returnable: subContent.aux?.returnable,
                                                chemicalData: chemicalData
                                              });
                                              
                                              if (chemicalData) {
                                                return renderChemicalStructures(chemicalData);
                                              }
                                              
                                              // Default rendering for non-chemical data
                                              return (
                                                <>
                                                  {subContent.content && (
                                                    <div className="sub-content-text">
                                                      {typeof subContent.content === 'string' ? 
                                                        subContent.content : 
                                                        <pre className="content-json">
                                                          {JSON.stringify(subContent.content, null, 2)}
                                                        </pre>
                                                      }
                                                    </div>
                                                  )}
                                                  {subContent.aux?.returnable && (
                                                    <div className="returnable-data">
                                                      <strong>Results:</strong>
                                                      <pre className="returnable-json">
                                                        {JSON.stringify(subContent.aux.returnable, null, 2)}
                                                      </pre>
                                                    </div>
                                                  )}
                                                </>
                                              );
                                            })()}
                                          </div>
                                        </div>
                                      ))}
                                    </div>
                                  </div>
                                )}
                              </div>
                              );
                            })}
                          </div>
                        ) : (
                          <div className="no-parsed-content">
                            <p>No reactions found in the session content.</p>
                            <p>The session content may not follow the expected format, or it may be empty.</p>
                          </div>
                        )}
                      </div>
                    ) : (
                      <pre className="content-json">
                        {JSON.stringify(session.content, null, 2)}
                      </pre>
                    )}
                  </div>
                ) : (
                  <div className="no-content">
                    <p>No content available yet.</p>
                    {session?.state === 'active' && (
                      <p className="content-hint">Content will appear here as the worker processes your query...</p>
                    )}
                  </div>
                )}
                
                {/* Output Files Section */}
                {(promptFile || csvFile) && (promptFile?.exists || csvFile?.exists) && (
                  <div className="output-files-section">
                    <div className="output-files-header">
                      <h3>Session Output Files</h3>
                      <span className="output-files-subtitle">Generated files from session processing</span>
                    </div>
                    <div className="output-files-grid">
                      {promptFile?.exists && (
                        <div className="output-file-item">
                          <div className="output-file-icon">
                            {React.createElement(MdTextSnippet as React.ComponentType<IconBaseProps>, { 
                              size: 24,
                              title: "Prompt Text File" 
                            })}
                          </div>
                          <div className="output-file-info">
                            <h4>Prompt File</h4>
                            <div className="output-file-details">
                              <span className="output-file-size">
                                {promptFile.content ? `${promptFile.content.length} characters` : 'Empty'}
                              </span>
                            </div>
                          </div>
                          <div className="output-file-actions">
                            {promptFile.content && (
                              <>
                                <button 
                                  className="output-file-btn view-btn"
                                  onClick={() => openPromptModal(promptFile.content!)}
                                  title="View prompt content"
                                >
                                  View
                                </button>
                                <button 
                                  className="output-file-btn download-btn"
                                  onClick={() => {
                                    // Create a blob and download
                                    const blob = new Blob([promptFile.content!], { type: 'text/plain' });
                                    const url = URL.createObjectURL(blob);
                                    const a = document.createElement('a');
                                    a.href = url;
                                    a.download = 'master_prompt.txt';
                                    document.body.appendChild(a);
                                    a.click();
                                    document.body.removeChild(a);
                                    URL.revokeObjectURL(url);
                                  }}
                                  title="Download prompt file"
                                >
                                  Download
                                </button>
                              </>
                            )}
                          </div>
                        </div>
                      )}
                      
                      {csvFile?.exists && (
                        <div className="output-file-item">
                          <div className="output-file-icon">
                            {React.createElement(MdTableChart as React.ComponentType<IconBaseProps>, { 
                              size: 24,
                              title: "CSV Data File" 
                            })}
                          </div>
                          <div className="output-file-info">
                            <h4>CSV Data</h4>
                            <div className="output-file-details">
                              <span className="output-file-size">
                                {csvFile.row_count} rows
                              </span>
                              {csvFile.headers && (
                                <span className="output-file-columns">
                                  {csvFile.headers.length} columns
                                </span>
                              )}
                            </div>
                          </div>
                          <div className="output-file-actions">
                            {csvFile.data && csvFile.headers && (
                              <>
                                <button 
                                  className="output-file-btn view-btn"
                                  onClick={() => openCsvModal(csvFile.data!, csvFile.headers!)}
                                  title="View CSV data table"
                                >
                                  View Table
                                </button>
                                <button 
                                  className="output-file-btn merge-btn"
                                  onClick={mergeCsvToProject}
                                  disabled={isMergingCsv}
                                  title="Merge session CSV with project master CSV"
                                >
                                  {isMergingCsv ? 'Merging...' : 'Merge to Project'}
                                </button>
                                <button 
                                  className="output-file-btn download-btn"
                                  onClick={() => {
                                    // Convert data back to CSV format
                                    const csvContent = [
                                      csvFile.headers!.join(','),
                                      ...csvFile.data!.map(row => 
                                        csvFile.headers!.map(header => 
                                          JSON.stringify(row[header] || '')
                                        ).join(',')
                                      )
                                    ].join('\n');
                                    
                                    const blob = new Blob([csvContent], { type: 'text/csv' });
                                    const url = URL.createObjectURL(blob);
                                    const a = document.createElement('a');
                                    a.href = url;
                                    a.download = 'data.csv';
                                    document.body.appendChild(a);
                                    a.click();
                                    document.body.removeChild(a);
                                    URL.revokeObjectURL(url);
                                  }}
                                  title="Download CSV file"
                                >
                                  Download
                                </button>
                              </>
                            )}
                          </div>
                        </div>
                      )}
                    </div>
                    
                    {loadingOutputFiles && (
                      <div className="output-files-loading">
                        <p>Loading output files...</p>
                      </div>
                    )}
                  </div>
                )}
                
                {/* Toggle All Button for Parsed View */}
                {viewMode === 'parsed' && parsedContent.length > 0 && (
                  <div className="toggle-all-container">
                    <button 
                      className="toggle-all-btn"
                      onClick={toggleAllReactions}
                    >
                      {expandedReactions.size === parsedContent.length ? 'Collapse All' : 'Expand All'}
                    </button>
                  </div>
                )}
              </div>
            </div>
          </>
        )}
      </div>
      
      {/* Individual Chemical Structure Fullscreen Modal */}
      {individualModal && (
        <div className="individual-modal-overlay" onClick={closeIndividualModal}>
          <div className="individual-modal-content" onClick={(e) => e.stopPropagation()}>
            <div className="individual-modal-header">
              <div className="individual-modal-title">
                <h2>#{individualModal.index + 1} - {(() => {
                  switch (individualModal.infoType) {
                    case 'smiles_update': return 'Composition SMILES';
                    case 'smarts_rxn_update': return 'Reaction SMARTS';
                    case 'smirks_rxn_update': return 'Reaction SMIRKS';
                    default: return individualModal.infoType.replace(/_/g, ' ').toUpperCase();
                  }
                })()}</h2>
              </div>
              <button 
                className="individual-close-btn"
                onClick={closeIndividualModal}
                title="Close fullscreen"
              >
                {React.createElement(MdClose as React.ComponentType<IconBaseProps>, { size: 24 })}
              </button>
            </div>
            <div className="individual-modal-body">
              <div className="individual-structure-container">
                <div className="individual-structure-input-section">
                  <label className="individual-structure-label">
                    {individualModal.type} Structure:
                  </label>
                  <textarea
                    className="individual-structure-textarea"
                    value={getCurrentStructureValue(individualModal.structure, individualModal.index, individualModal.infoType)}
                    onChange={(e) => handleStructureChange(individualModal.structureKey, e.target.value)}
                    placeholder={`Enter ${individualModal.type} structure`}
                    rows={6}
                    autoFocus
                  />
                </div>
                <div className="individual-structure-display">
                  <ChemicalRenderer
                    data={getCurrentStructureValue(individualModal.structure, individualModal.index, individualModal.infoType)}
                    inputType={individualModal.type}
                    className="fullscreen-single"
                    onRenderError={(error) => console.warn(`Failed to render ${individualModal.type}:`, error)}
                  />
                </div>
              </div>
            </div>
          </div>
        </div>
      )}
      
      {/* Prompt File Modal */}
      {promptModal.isOpen && (
        <div className="modal-overlay" onClick={closePromptModal}>
          <div className="modal-content prompt-modal" onClick={(e) => e.stopPropagation()}>
            <div className="modal-header">
              <div className="modal-title">
                <h2>Prompt File Content</h2>
              </div>
              <button 
                className="modal-close-btn"
                onClick={closePromptModal}
                title="Close"
              >
                {React.createElement(MdClose as React.ComponentType<IconBaseProps>, { size: 24 })}
              </button>
            </div>
            <div className="modal-body">
              <div className="prompt-content-container">
                <textarea
                  className="prompt-content-textarea"
                  value={promptModal.content}
                  readOnly
                  rows={20}
                />
              </div>
            </div>
            <div className="modal-footer">
              <button
                className="modal-btn download-btn"
                onClick={() => {
                  const blob = new Blob([promptModal.content], { type: 'text/plain' });
                  const url = URL.createObjectURL(blob);
                  const a = document.createElement('a');
                  a.href = url;
                  a.download = 'master_prompt.txt';
                  document.body.appendChild(a);
                  a.click();
                  document.body.removeChild(a);
                  URL.revokeObjectURL(url);
                }}
              >
                Download
              </button>
            </div>
          </div>
        </div>
      )}
      
      {/* CSV Table Modal */}
      {csvModal.isOpen && (
        <div className="modal-overlay" onClick={closeCsvModal}>
          <div className="modal-content csv-modal" onClick={(e) => e.stopPropagation()}>
            <div className="modal-header">
              <div className="modal-title">
                <h2>CSV Data Table</h2>
                <span className="csv-info">
                  {csvModal.data.length} rows × {csvModal.headers.length} columns
                </span>
              </div>
              <div className="modal-actions">
                <button
                  className={`modal-btn edit-toggle-btn ${csvModal.isEditing ? 'active' : ''}`}
                  onClick={toggleCsvEditing}
                  disabled={isSavingCsv}
                >
                  {csvModal.isEditing ? 'View Mode' : 'Edit Mode'}
                </button>
                <button 
                  className="modal-close-btn"
                  onClick={closeCsvModal}
                  title="Close"
                >
                  {React.createElement(MdClose as React.ComponentType<IconBaseProps>, { size: 24 })}
                </button>
              </div>
            </div>
            <div className="modal-body csv-modal-body">
              <div className="csv-table-container">
                <table className="csv-table">
                  <thead>
                    <tr>
                      {csvModal.headers.map((header, index) => (
                        <th key={index} className="csv-header">
                          {header}
                        </th>
                      ))}
                      {csvModal.isEditing && (
                        <th className="csv-header csv-actions-header">
                          Actions
                        </th>
                      )}
                    </tr>
                  </thead>
                  <tbody>
                    {(csvModal.isEditing ? editedCsvData : csvModal.data).map((row, rowIndex) => (
                      <tr key={rowIndex} className="csv-row">
                        {csvModal.headers.map((header, colIndex) => (
                          <td key={colIndex} className="csv-cell">
                            {csvModal.isEditing ? (
                              <input
                                type="text"
                                className="csv-cell-input"
                                value={editedCsvData[rowIndex]?.[header] || ''}
                                onChange={(e) => handleCsvCellChange(rowIndex, header, e.target.value)}
                              />
                            ) : (
                              <span className="csv-cell-value">
                                {row[header] || ''}
                              </span>
                            )}
                          </td>
                        ))}
                        {csvModal.isEditing && (
                          <td className="csv-cell csv-actions-cell">
                            <button
                              className="csv-delete-btn"
                              onClick={() => deleteCsvRow(rowIndex)}
                              title="Delete this row"
                              disabled={isSavingCsv}
                            >
                              {React.createElement(MdDelete as React.ComponentType<IconBaseProps>, { size: 16 })}
                            </button>
                          </td>
                        )}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
            <div className="modal-footer">
              {csvModal.isEditing && (
                <div className="edit-actions">
                  <button
                    className="modal-btn save-btn"
                    onClick={saveCsvChanges}
                    disabled={isSavingCsv}
                  >
                    {isSavingCsv ? 'Saving...' : 'Save Changes'}
                  </button>
                  <button
                    className="modal-btn reset-btn"
                    onClick={resetCsvChanges}
                    disabled={isSavingCsv}
                  >
                    Reset
                  </button>
                </div>
              )}
              <button
                className="modal-btn download-btn"
                onClick={() => {
                  const dataToDownload = csvModal.isEditing ? editedCsvData : csvModal.data;
                  const csvContent = [
                    csvModal.headers.join(','),
                    ...dataToDownload.map(row => 
                      csvModal.headers.map(header => 
                        JSON.stringify(row[header] || '')
                      ).join(',')
                    )
                  ].join('\n');
                  
                  const blob = new Blob([csvContent], { type: 'text/csv' });
                  const url = URL.createObjectURL(blob);
                  const a = document.createElement('a');
                  a.href = url;
                  a.download = 'data.csv';
                  document.body.appendChild(a);
                  a.click();
                  document.body.removeChild(a);
                  URL.revokeObjectURL(url);
                }}
              >
                Download CSV
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default SessionDetails;
