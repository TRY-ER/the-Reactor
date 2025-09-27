import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import './index.css';
import { Dashboard, ProjectDetails, SessionDetails } from './components';

function App() {
  return (
    <div className="App">
      <Router>
        <Routes>
          <Route path="/" element={<Dashboard />} />
          <Route path="/project/:id" element={<ProjectDetails />} />
          <Route path="/project/:projectId/session/:sessionId" element={<SessionDetails />} />
        </Routes>
      </Router>
    </div>
  );
}

export default App;
